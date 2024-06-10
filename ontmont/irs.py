import random
from operator import xor

import pandas as pd

from .utils import reverse_complement, shuffle_seq


def calc_distance_score(dist1, dist2, dist_cutoff1, dist_cutoff2):
    dist_score = ( 2 * (dist1 <= dist_cutoff1) + (dist1 <= dist_cutoff2) # at least one should be <= cutoff1
                 + 2 * (dist2 <= dist_cutoff1) + (dist2 <= dist_cutoff2) )
    return dist_score

def align_two_sequences(seq1, seq2, sw, rc=True):
    len1, len2 = len(seq1), len(seq2)
    seq2_rc = reverse_complement(seq2)
    if rc:
        aln = sw.align(seq1, seq2_rc)
    else:
        aln = sw.align(seq1, seq2)
    aln.seq1 = seq1
    aln.seq2 = seq2
    down_start = (len2 - aln.q_end) if rc else aln.q_pos
    down_end = (len2 - aln.q_pos) if rc else aln.q_end
    up_start = aln.r_pos
    up_end = aln.r_end
    up = (up_start, up_end)
    down = (down_start, down_end)
    aln.up_coords = up
    aln.down_coords = down
    return aln


def overlaps(seg1, seg2):
    start1, end1 = seg1
    start2, end2 = seg2
    if (start1 < end2) and (start2 < end1):
        return True
    return False

def calc_pval_onesided_ir(seq, sw, direction, src_score, dist_cutoff, n_iter=1000, random_seed=42):
    random.seed(random_seed)
    scores = []
    hairpin_flags = []
    length = len(seq)
    for i in range(n_iter):
        hairpin_flag = 0
        _seq = shuffle_seq(seq)
        _aln = get_onesided_ir(_seq, sw)
        up_start = length - _aln.q_end
        up_end = length - _aln.q_pos
        down_start = _aln.r_pos
        down_end = _aln.r_end
        up = (up_start, up_end)
        down = (down_start, down_end)
        if not overlaps(up, down):
            dist = (up_start) if direction=='down' else (length - down_end)
            score = _aln.score
            scores.append(score)
            is_near_breakpoint = (dist <= dist_cutoff)
            if score >= src_score and is_near_breakpoint:
                hairpin_flag = 1
        hairpin_flags.append(hairpin_flag)
    scores = pd.Series(scores)
    hairpin_flags = pd.Series(hairpin_flags)
    pval = hairpin_flags.mean()
    return pval

def get_best_onesided_ir(seq, direction, sw, dist_cutoff=1, margins=[15, 30, 60]):
    max_margin = max(margins)
    assert len(seq) >= max_margin
    _seq = seq
    best_score = 0
    best_aln = None
    best_pvalue = 1
    for margin in margins:
        seq = _seq[:margin] if direction=='down' else _seq[-margin:]
        length = len(seq)
        aln = get_onesided_ir(seq, sw)
        up_start = length - aln.q_end
        up_end = length - aln.q_pos
        down_start = aln.r_pos
        down_end = aln.r_end
        up = (up_start, up_end)
        down = (down_start, down_end)
        if overlaps(up, down):
            continue
        aln.up_coords = up
        aln.down_coords = down
        aln.seq = seq
        dist = (up_start) if direction=='down' else (length - down_end)
        # print(margin, dist, (aln.r_pos), (aln.q_pos), best_score, dist_cutoff)

        if aln.score > best_score:
            if dist <= dist_cutoff:
                best_score = aln.score
                best_aln = aln

    if best_aln:
        pvalue = calc_pval_onesided_ir(best_aln.seq, sw, direction, 
            best_score, dist_cutoff, n_iter=100, random_seed=42)
        best_aln.pvalue = pvalue
    return best_aln
        
def get_onesided_ir(seq, sw):
    length = len(seq)
    seq_rc = reverse_complement(seq)
    aln = sw.align(seq, seq_rc)
    up_start = length - aln.q_end
    up_end = length - aln.q_pos
    down_start = aln.r_pos
    down_end = aln.r_end
    up = (up_start, up_end)
    aln.up_coords = up
    down = (down_start, down_end)
    aln.down_coords = down
    return aln
    

def calc_pval_bilateral_ir(seq1, seq2, sw, src_score, dist_cutoff1, dist_cutoff2, n_iter=1000, random_seed=42):
    dist_score_sum = 4 # hard coded
    random.seed(random_seed)
    hairpin_flags = []

    for i in range(n_iter):
        hairpin_flag = 0
        _seq1 = shuffle_seq(seq1)
        _seq2 = shuffle_seq(seq2)
        _aln = align_two_sequences(_seq1, _seq2, sw, rc=True)
        _aln.dist1 = len(seq1) - _aln.up_coords[1] # up end
        _aln.dist2 = _aln.down_coords[0] # down start
    
        if _aln.score >= src_score:
            _aln.dist_score = calc_distance_score(_aln.dist1, _aln.dist2, dist_cutoff1, dist_cutoff2)
            if _aln.dist_score >= dist_score_sum:
                hairpin_flag = 1
        hairpin_flags.append(hairpin_flag)
    hairpin_flags = pd.Series(hairpin_flags)
    
    pval = hairpin_flags.mean()
    return pval


def get_best_ir_within_breakpoints(seq1, seq2, sw, dist_cutoff1=1, dist_cutoff2=3, margins=[15, 30, 60]):
    dist_score_sum = 4
    max_margin = max(margins)
    best_score = 0
    best_aln = None
    assert len(seq1) >= max_margin or len(seq2) >= max_margin
    _seq1, _seq2 = seq1, seq2
    for margin in margins:
        seq1 = _seq1[-margin:]
        seq2 = _seq2[:margin]
        len1 = len(seq1)
        aln = align_two_sequences(seq1, seq2, sw)
        aln.dist1 = len1 - aln.up_coords[1] # up end
        aln.dist2 = aln.down_coords[0] # down start

        if aln.score > best_score:
            aln.dist_score = calc_distance_score(aln.dist1, aln.dist2, dist_cutoff1, dist_cutoff2)
            if aln.dist_score >= dist_score_sum:
                best_score = aln.score
                best_aln = aln

    if best_aln:
        pvalue = calc_pval_bilateral_ir(best_aln.seq1, best_aln.seq2, sw, best_score, 
            dist_cutoff1, dist_cutoff2, n_iter=100, random_seed=42)
        best_aln.pvalue = pvalue

    return best_aln


def calc_pval_segmental_ir(seq1, seq2, ori1, ori2, sw, src_score, dist_cutoff1, dist_cutoff2, n_iter=1000, random_seed=42):
    dist_score_sum = 4 # hard coded
    random.seed(random_seed)
    hairpin_flags = []
    len1, len2 = len(seq1), len(seq2)
    seq = seq1 + seq2

    for i in range(n_iter):
        hairpin_flag = 0
        _seq1 = shuffle_seq(seq1)
        _seq2 = shuffle_seq(seq2)
        _aln = align_two_sequences(_seq1, _seq2, sw, rc=True)
        _aln.dist1 = len1 - _aln.up_coords[1] if ori1 == '+' else _aln.up_coords[0] # seq1
        _aln.dist2 = len2 - _aln.down_coords[1] if ori2 == '+' else _aln.down_coords[0] # seq2

    
        if _aln.score >= src_score:
            _aln.dist_score = calc_distance_score(_aln.dist1, _aln.dist2, dist_cutoff1, dist_cutoff2)
            if _aln.dist_score >= dist_score_sum:
                hairpin_flag = 1
        hairpin_flags.append(hairpin_flag)
    hairpin_flags = pd.Series(hairpin_flags)
    
    pval = hairpin_flags.mean()
    return pval


def get_best_ir_within_segment(pair, sw, genome, dist_cutoff1=1, dist_cutoff2=3, margins=[15, 30, 60]):
    dist_score_sum = 4
    best_score = 0
    best_aln = None
    for margin in margins:
        pair.brk1.get_breakpoint_seqs(margin, genome)
        pair.brk2.get_breakpoint_seqs(margin, genome)
        seq1 = pair.brk1.seq_rearranged # seq that remains
        seq2 = pair.brk2.seq_rearranged # seq that remains
        len1, len2 = len(seq1), len(seq2)
        aln = align_two_sequences(seq1, seq2, sw)
        aln.flag_segment = False
        aln.dist1 = len1 - aln.up_coords[1] if pair.brk1.ori == '+' else aln.up_coords[0] # seq1
        aln.dist2 = len2 - aln.down_coords[1] if pair.brk2.ori == '+' else aln.down_coords[0] # seq2

        if aln.score > best_score:
            aln.dist_score = calc_distance_score(aln.dist1, aln.dist2, dist_cutoff1, dist_cutoff2)
            if aln.dist_score >= dist_score_sum:
                best_score = aln.score
                aln.flag_segment = True
                best_aln = aln

    if best_aln:
        pvalue = calc_pval_segmental_ir(best_aln.seq1, best_aln.seq2, 
            pair.brk1.ori, pair.brk2.ori, 
            sw, best_score, 
            dist_cutoff1, dist_cutoff2, n_iter=100, random_seed=42)
        best_aln.pvalue = pvalue

    pair.aln_segment = best_aln
    return pair


def is_holliday_junction(aln, negative1, negative2, 
                         dist_cutoff1=1, dist_cutoff2=3, dist_sum_cutoff=4, score_cutoff=4):
    len1, len2 = len(aln.seq1), len(aln.seq2)
    aln.dist1 = (aln.r_pos) if negative1 else (len1 - aln.r_end)
    aln.dist2 = (aln.q_pos) if negative2 else (len2 - aln.q_end)
    aln.dist_score = calc_distance_score(aln.dist1, aln.dist2, dist_cutoff1, dist_cutoff2)
    aln.flag_holliday = False
    if aln.dist_score >= dist_sum_cutoff and aln.score >= score_cutoff:
        aln.flag_holliday = True
    return aln


def get_breakpoint_pair_seq_data(pair):
    seqs = {}
    up1, dn1 = pair.brk1.upstream, pair.brk1.downstream
    up2, dn2 = pair.brk2.upstream, pair.brk2.downstream
    seqs = {'u1':up1, 'u2':up2, 'd1':dn1, 'd2':dn2}
    seqs['u2r'] = seqs['u2']
    seqs['d2r'] = seqs['d2']
    return seqs

def calc_pval_holliday(seq1, seq2, negative1, negative2, src_score, sw, 
                       dist_cutoff1=1, dist_cutoff2=3, n_iter=1000, random_seed=42):
    dist_sum_cutoff = 4 # hard coded
    random.seed(random_seed)
    flags = []

    for i in range(n_iter):
        flag = 0
        _seq1 = shuffle_seq(seq1)
        _seq2 = shuffle_seq(seq2)
        _aln = align_two_sequences(_seq1, _seq2, sw, rc=False)
        _aln = is_holliday_junction(_aln, negative1, negative2, 
            dist_cutoff1, dist_cutoff2, 
            dist_sum_cutoff=dist_sum_cutoff, score_cutoff=src_score)
        if _aln.score > src_score and _aln.flag_holliday:
            flag = 1
        
        flags.append(flag)
    flags = pd.Series(flags)
    pval = flags.mean()
    return pval

def get_best_holliday_junctions(pair, sw, genome, score_cutoff=4, dist_cutoff1=2, dist_cutoff2=5, margins=[15, 30, 60]):
    dist_sum_cutoff = 4 # HARD CODED CUTOFF
    key_pairs = [('u1', 'u2'), ('d1', 'd2'), ('u1', 'd2r'), ('d1', 'u2r')]
    alns = {}
    # pvals = {}
    for src, dst in key_pairs:
        src_negative = xor(src.startswith('d'), src.endswith('r'))
        dst_negative = xor(dst.startswith('d'), dst.endswith('r'))
        dst_reverse = dst[-1] == 'r'

        best_score = 0
        best_aln = None
        for margin in margins:
            pair.brk1.get_breakpoint_seqs(margin, genome)
            pair.brk2.get_breakpoint_seqs(margin, genome)
            seqs = get_breakpoint_pair_seq_data(pair)

            key = f'{src}_{dst}'
            aln = align_two_sequences(seqs[src], seqs[dst], sw, rc=dst_reverse)
            src_negative = xor(src.startswith('d'), src.endswith('r'))
            dst_negative = xor(dst.startswith('d'), dst.endswith('r'))
            aln = is_holliday_junction(aln, src_negative, dst_negative, 
                                       dist_cutoff1, dist_cutoff2, dist_sum_cutoff, score_cutoff)
            if aln.flag_holliday and aln.score > best_score:
                best_score = aln.score
                best_aln = aln
        if best_aln:
            pvalue = calc_pval_holliday(seqs[src], seqs[dst], src_negative, dst_negative, 
                best_aln.score, sw, dist_cutoff1, dist_cutoff2, 
                n_iter=100, random_seed=42)
            best_aln.pvalue = pvalue

        alns[key] = best_aln

    pair.alns = alns
    return pair

