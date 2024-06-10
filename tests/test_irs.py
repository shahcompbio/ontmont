import yaml
import swalign
import pyfaidx 
import pytest

from ontmont.datatypes import (Breakpoint, BreakpointPair)
from ontmont.utils import reverse_complement
from ontmont.irs import (calc_distance_score, get_best_onesided_ir, align_two_sequences,
    get_onesided_ir, get_best_ir_within_breakpoints, get_best_ir_within_segment,
    calc_pval_onesided_ir, calc_pval_bilateral_ir, calc_pval_segmental_ir, calc_pval_holliday,
    overlaps, is_holliday_junction, get_best_holliday_junctions, get_breakpoint_pair_seq_data)

@pytest.fixture
def config():
    CONFIG = """
    palindrome:
        match: 1
        mismatch: -2
        gap_penalty: -3
        gap_extension_penalty: -1
        cutoff1: 6
        cutoff2: 8
        bp_near_breakpoint: 5
        random_seed: 42
        n_iter: 1000
        margin: 30

    holliday:
        match: 1
        mismatch: -2
        gap_penalty: -2
        gap_extension_penalty: -1
        distance_cutoff1: 1
        distance_cutoff2: 3
        score_cutoff: 5
        random_seed: 42
        n_iter: 1000
        margin: 15

    reference:
        '1': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_GFPVector.fa
        '2': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_PGBD5Vector.fa
        '3': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_PiggyBacVector.fa
        '5': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_DelPBEF1_PGBD5Vector.fa
    """
    return yaml.safe_load(CONFIG)


@pytest.fixture
def sw_palindrome(config):
    scoring = swalign.NucleotideScoringMatrix(float(config['palindrome']['match']), float(config['palindrome']['mismatch']))
    sw = swalign.LocalAlignment(scoring, globalalign=False,
        gap_penalty=float(config['palindrome']['gap_penalty']), 
        gap_extension_penalty=float(config['palindrome']['gap_extension_penalty']))
    return sw

@pytest.fixture
def sw_holliday(config):
    scoring = swalign.NucleotideScoringMatrix(float(config['holliday']['match']), float(config['holliday']['mismatch']))
    sw = swalign.LocalAlignment(scoring, globalalign=False,
        gap_penalty=float(config['holliday']['gap_penalty']), 
        gap_extension_penalty=float(config['holliday']['gap_extension_penalty']))
    return sw

# TESTS

def test_reverse_complement():
    seq_revs = {'ACGGT':'ACCGT', 'GNTTCA': 'TGAA.C'}
    for seq, rev in seq_revs.items():
        assert reverse_complement(seq) == rev, (seq, rev)

def test_calc_distance_score():
    assert calc_distance_score(dist1=0, dist2=0, dist_cutoff1=1, dist_cutoff2=3) == 6
    assert calc_distance_score(dist1=1, dist2=3, dist_cutoff1=1, dist_cutoff2=3) == 4
    assert calc_distance_score(dist1=2, dist2=3, dist_cutoff1=1, dist_cutoff2=3) == 2
    assert calc_distance_score(dist1=2, dist2=4, dist_cutoff1=1, dist_cutoff2=3) == 1
    assert calc_distance_score(dist1=4, dist2=4, dist_cutoff1=1, dist_cutoff2=3) == 0


def test_align_two_sequences(config, sw_palindrome):
    genome = pyfaidx.Fasta(config['reference']['3'])

    brk = Breakpoint('PBEF1NeoTransposon', 1478, '+')
    brk.get_breakpoint_seqs(margin=50, genome=genome)
    seq1 = brk.upstream
    seq2 = brk.downstream
    sw = sw_palindrome
    aln = align_two_sequences(seq1, seq2, sw, rc=True)
    up_start, up_end = aln.up_coords
    down_start, down_end = aln.down_coords
    assert aln.seq1[up_start:up_end] == 'ATGCA', aln.seq1[up_start:up_end]
    assert aln.seq2[down_start:down_end] == reverse_complement('ATGCA'), aln.dump()

    seq1 = 'TTAAAAAAAT'
    seq2 = 'CAAAAAAACC'
    aln = align_two_sequences(seq1, seq2, sw, rc=False) # not rc
    assert aln.up_coords == (2, 2+7), aln.up_coords
    assert aln.down_coords == (1, 1+7), aln.down_coords
    up_start, up_end = aln.up_coords
    down_start, down_end = aln.down_coords
    assert aln.seq1[up_start:up_end] == 'AAAAAAA', aln.seq1[up_start:up_end]
    assert aln.seq2[down_start:down_end] == 'AAAAAAA', aln.seq2[down_start:down_end]


def test_get_best_onesided_ir(sw_palindrome):
    seq = 'TTTTTTTTTTTTAGCACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGCTAA'
    sw = sw_palindrome
    direction = 'up' # sequence is upstream of a breakpoint

    aln = get_best_onesided_ir(seq, direction, sw, margins=[15, 30, 60])
    assert aln.score == 7, aln.dump()
    assert aln.seq[aln.up_coords[0]:aln.up_coords[1]] == 'TTAGCAC', aln.dump()
    assert aln.seq[aln.down_coords[0]:aln.down_coords[1]] == 'GTGCTAA', aln.dump()

    aln = get_best_onesided_ir(seq, direction, sw, margins=[15, 30])
    assert aln.score == 2, aln.dump()
    assert aln.seq[aln.up_coords[0]:aln.up_coords[1]] == 'TT', aln.dump()
    assert aln.seq[aln.down_coords[0]:aln.down_coords[1]] == 'AA', aln.dump()

    seq = 'GAGCTTGTTGTACAG'
    direction = 'up'
    aln = get_best_onesided_ir(seq, direction, sw, margins=[15]) # , ((up_start, up_end), (down_start, down_end))
    assert aln == None, aln

def test_get_onesided_ir(sw_palindrome):
    sw = sw_palindrome
    seq = 'TTAGCACTTTTTTTTTTTTTGTGCTAAAAA'
    aln = get_onesided_ir(seq, sw) # , ((up_start, up_end), (down_start, down_end))
    up_start, up_end = aln.up_coords
    down_start, down_end = aln.down_coords
    assert up_start == 0, aln.dump()
    assert up_end == 7, aln.dump()
    assert down_start == 20, aln.dump()
    assert down_end == 27, aln.dump()
    assert seq[up_start:up_end] == reverse_complement(seq[down_start:down_end]) == 'TTAGCAC', (
        seq[up_start:up_end], seq[down_start:down_end]
    )
    assert len(seq) - down_end == 3, aln.dump() 

def test_calc_pval_onesided_ir_1(sw_palindrome):
    sw = sw_palindrome
    seq = 'GGGGAAAAATTTTTGGGGGGGAAAAATTTTT'
    direction = 'up'
    aln = get_onesided_ir(seq, sw)
    src_score = aln.score
    pval = calc_pval_onesided_ir(seq, sw, direction, src_score, dist_cutoff=5, n_iter=100, random_seed=42)
    assert src_score == 10.0, src_score
    assert pval == 0.0, pval

def test_calc_pval_onesided_ir_2(sw_palindrome):
    sw = sw_palindrome
    seq = 'CCCTTGGGGAAGGGTGGGGTGACGCGGGGG' # downstream
    direction = 'down'
    aln = get_onesided_ir(seq, sw)
    src_score = aln.score
    pval = calc_pval_onesided_ir(seq, sw, direction, src_score, dist_cutoff=5, n_iter=100, random_seed=42)
    assert src_score == 5.0, src_score
    assert pval == 0.01, pval

def test_calc_pval_onesided_ir_3(sw_palindrome):
    sw = sw_palindrome
    seq = 'GGGGGCGCAGTGGGGTGGGAAGGGGTTCCC' # upstream, reverse of seq from _2
    direction = 'up'
    aln = get_onesided_ir(seq, sw)
    src_score = aln.score
    pval = calc_pval_onesided_ir(seq, sw, direction, src_score, dist_cutoff=5, n_iter=100, random_seed=42)
    assert src_score == 5.0, src_score
    assert pval == 0.01, pval

def test_overlaps():
    seg1 = (8, 14)
    seg2 = (8, 14)
    assert overlaps(seg1, seg2), (seg1, seg2)
    seg1 = (8, 14)
    seg2 = (13, 15)
    assert overlaps(seg1, seg2), (seg1, seg2)

def test_get_best_ir_within_breakpoints(sw_palindrome):
    seq1 = 'TTTTTTTTTTTTTCGCGAGCCGTCGTTGCGCTTGTGCTAA'
    seq2 = 'TTAGCACTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAA'
    sw = sw_palindrome

    aln = get_best_ir_within_breakpoints(seq1, seq2, sw, dist_cutoff1=1, dist_cutoff2=3, margins=[15, 40])
    assert aln.score == 7, aln.dump()
    assert aln.seq1[aln.up_coords[0]:aln.up_coords[1]] == 'GTGCTAA', aln.dump()
    assert aln.seq2[aln.down_coords[0]:aln.down_coords[1]] == 'TTAGCAC', aln.dump()

def test_calc_pval_bilateral_ir_1(sw_palindrome):
    sw = sw_palindrome
    seq1 = 'AAAAAAAAAA'
    seq2 = 'TTTTTTTTTT'
    aln = align_two_sequences(seq1, seq2, sw, rc=True)
    src_score = aln.score
    pval = calc_pval_bilateral_ir(seq1, seq2, sw, src_score, dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 1, pval

def test_calc_pval_bilateral_ir_2(sw_palindrome):
    sw = sw_palindrome
    seq1 = 'AAAAAAAAAA'
    seq2 = 'TTAAATTTTT'
    aln = align_two_sequences(seq1, seq2, sw, rc=True)
    src_score = aln.score
    pval = calc_pval_bilateral_ir(seq1, seq2, sw, src_score, dist_cutoff1=100, dist_cutoff2=10, n_iter=10, random_seed=42)
    assert pval == 0.3, pval

def test_get_best_ir_within_segment(config, sw_palindrome):
    genome = pyfaidx.Fasta(config['reference']['3'])
    sw = sw_palindrome
    margin = 50
    chrom1, pos1, ori1 = 'PBEF1NeoTransposon', 1478, '-'
    chrom2, pos2, ori2 = 'PBEF1NeoTransposon', 4996, '+'
    brk1 = Breakpoint(chrom1, pos1, ori1)
    brk2 = Breakpoint(chrom2, pos2, ori2)
    pair = BreakpointPair(brk1, brk2)
    pair.brk1.get_breakpoint_seqs(margin, genome)
    pair.brk2.get_breakpoint_seqs(margin, genome)
    
    pair = get_best_ir_within_segment(pair, sw, genome, margins=[margin])
    assert pair.aln_segment.flag_segment, pair.aln_segment.dump()
    assert pair.aln_segment.score == 17, pair.aln_segment.dump()
    start1, end1 = pair.aln_segment.up_coords
    assert pair.brk1.seq_rearranged[start1:end1] == 'TTAACCCTAGAAAGATAGTC', pair.brk1.seq_rearranged[start1:end1]
    start2, end2 = pair.aln_segment.down_coords
    assert pair.brk2.seq_rearranged[start2:end2] == 'GATTATCTTTCTAGGGTTAA', pair.brk2.seq_rearranged[start2:end2]

def test_calc_pval_segmental_ir_1(sw_palindrome):
    sw = sw_palindrome
    seq1 = 'AAAAAAAAAA'
    seq2 = 'TTTTTTTTTT'
    ori1 = '-'
    ori2 = '+'
    aln = align_two_sequences(seq1, seq2, sw, rc=True)
    aln.dump()
    src_score = aln.score
    pval = calc_pval_segmental_ir(seq1, seq2, ori1, ori2, sw, src_score, 
                                  dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 1, pval # shuffling sequence doesn't affect aln

def test_calc_pval_segmental_ir_2(sw_palindrome):
    sw = sw_palindrome
    seq1 = 'AAAAGGCGG'
    seq2 = 'TTAAATTTT'
    ori1 = '-'
    ori2 = '+'
    aln = align_two_sequences(seq1, seq2, sw, rc=True)
    src_score = aln.score
    pval = calc_pval_segmental_ir(seq1, seq2, ori1, ori2, sw, src_score, 
                                  dist_cutoff1=1000, dist_cutoff2=10, n_iter=10, random_seed=42)
    assert pval == 0.1, pval

def test_is_holliday_junction(sw_palindrome):
    score_cutoff = 4
    dist_sum_cutoff = 4 # HARD CODED
    sw = sw_palindrome
    seq1 = 'AACCGGTTACGTACGT'
    seq2 = 'TTTTACGTAACCGGTA'
    negative1, negative2 = False, False
    aln = align_two_sequences(seq1, seq2, sw, rc=False)
    aln = is_holliday_junction(aln, negative1, negative2, 
        dist_cutoff1=1, dist_cutoff2=9, dist_sum_cutoff=dist_sum_cutoff, score_cutoff=score_cutoff)
    assert aln.flag_holliday, aln.dump()

def test_calc_pval_holliday_u1u2(sw_holliday):
    sw = sw_holliday
    seq1 = 'TTACGC'+'CAAGCGCGC'
    seq2 = 'GCACAG'+'CAAGCGCGC'
    negative1 = False # downstream ^ rc?
    negative2 = False # downstream ^ rc?
    aln = align_two_sequences(seq1, seq2, sw, rc=False) # no rc
    assert aln.score == 9, aln.dump()
    src_score = aln.score
    pval = calc_pval_holliday(seq1, seq2, negative1, negative2, src_score, sw,
                              dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 0.0, pval

def test_calc_pval_holliday_d1d2(sw_holliday):
    sw = sw_holliday
    seq1 = 'CAAGC'+'TTACGC'
    seq2 = 'CAAGC'+'GCACAG'
    negative1 = True # downstream ^ rc?
    negative2 = True # downstream ^ rc?
    aln = align_two_sequences(seq1, seq2, sw, rc=False) # no rc
    assert aln.score == 5, aln.dump()
    src_score = aln.score
    pval = calc_pval_holliday(seq1, seq2, negative1, negative2, src_score, sw,
                              dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 0.02, pval

def test_calc_pval_holliday_u1d2r(sw_holliday):
    sw = sw_holliday
    seq1 = 'TTACGC'+'CAAG'
    seq2 = 'CAAG'+'GCACAG'
    negative1 = False # downstream ^ rc?
    negative2 = True # downstream ^ rc?
    aln = align_two_sequences(seq1, seq2, sw, rc=False) # no rc
    aln = is_holliday_junction(aln, negative1, negative2, 
                               dist_cutoff1=1, dist_cutoff2=3, 
                               dist_sum_cutoff=4, score_cutoff=4)
    assert aln.flag_holliday, aln.dump()
    assert aln.score == 4, aln.dump()
    src_score = aln.score
    pval = calc_pval_holliday(seq1, seq2, negative1, negative2, src_score, sw,
                              dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 0.01, pval

def test_calc_pval_holliday_d1u2r(sw_holliday):
    sw = sw_holliday
    seq1 = 'CAAG'+'GCACAG'
    seq2 = 'TTACGC'+'CAAG'
    negative1 = True # downstream ^ rc?
    negative2 = False # downstream ^ rc?
    aln = align_two_sequences(seq1, seq2, sw, rc=False) # no rc
    aln = is_holliday_junction(aln, negative1, negative2, 
                               dist_cutoff1=1, dist_cutoff2=3, 
                               dist_sum_cutoff=4, score_cutoff=4)
    assert aln.flag_holliday, aln.dump()
    assert aln.score == 4, aln.dump()
    src_score = aln.score
    pval = calc_pval_holliday(seq1, seq2, negative1, negative2, src_score, sw,
                              dist_cutoff1=1, dist_cutoff2=3, n_iter=100, random_seed=42)
    assert pval == 0.01, pval

def test_get_breakpoint_pair_seq_data():
    brk1 = Breakpoint('chr1', 1000, '+')
    brk2 = Breakpoint('chrX', 2000, '-')
    brk1.upstream = 'AAAAA' # up1
    brk1.downstream = 'CCCCG' # dn1
    brk2.upstream = 'GCGGG' # up2
    brk2.downstream = 'ATTTT' # dn2
    pair = BreakpointPair(brk1, brk2)
    seqs = get_breakpoint_pair_seq_data(pair)
    assert seqs['u1'] == 'AAAAA', seqs['u1']
    assert seqs['u2'] == 'GCGGG', seqs['u2']
    assert seqs['d1'] == 'CCCCG', seqs['d1']
    assert seqs['d2'] == 'ATTTT', seqs['d2']
    assert seqs['u2r'] == 'GCGGG', seqs['u2r']
    assert seqs['d2r'] == 'ATTTT', seqs['d2r']


def test_get_best_holliday_junction_1(config, sw_holliday): 
    genome = pyfaidx.Fasta(config['reference']['2'])
    sw = sw_holliday

    chrom1, pos1, ori1 = 'PGBD5Vector', 1012, '+'
    chrom2, pos2, ori2 = 'PGBD5Vector', 1980, '+'
    brk1 = Breakpoint(chrom1, pos1, ori1)
    brk2 = Breakpoint(chrom2, pos2, ori2)
    pair = BreakpointPair(brk1, brk2)
    pair = get_best_holliday_junctions(pair, sw, genome)
    alns = pair.alns
    assert alns['u1_u2'].score == 6, alns['u1_u2'].dump()
    start1, end1 = alns['u1_u2'].up_coords
    assert (start1, end1) == (6, 15), (start1, end1)
    assert alns['u1_u2'].seq1[start1:end1] == 'CAAGCGCGC', alns['u1_u2'].seq1
    start2, end2 = alns['u1_u2'].down_coords
    assert (start2, end2) == (6, 14), (start2, end2)
    assert alns['u1_u2'].seq2[start2:end2] == 'CAAGCG-GC'.replace('-', ''), alns['u1_u2'].seq1

def test_get_best_holliday_junction_2(config, sw_holliday): 
    genome = pyfaidx.Fasta(config['reference']['2'])
    sw = sw_holliday

    chrom1, pos1, ori1 = 'PGBD5Vector', 1285, '-'
    chrom2, pos2, ori2 = 'PGBD5Vector', 3421, '-'
    brk1 = Breakpoint(chrom1, pos1, ori1)
    brk2 = Breakpoint(chrom2, pos2, ori2)
    pair = BreakpointPair(brk1, brk2)
    pair = get_best_holliday_junctions(pair, sw, genome)
    alns = pair.alns
    
    assert alns['u1_u2'] == None, alns['u1_u2']
    assert alns['d1_d2'] == None, alns['d1_d2']
    assert alns['u1_d2r'] == None, alns['u1_d2r']
    assert alns['d1_u2r'].score == 5, alns['d1_u2r'].dump()

    start1, end1 = alns['d1_u2r'].up_coords
    assert (start1, end1) == (0, 5), (start1, end1)
    assert pair.brk1.downstream[start1:end1] == 'GGGTC', alns['d1_u2r'].dump() # pair.brk1.downstream #pair.brk1.upstream[start1:end1]
    start2, end2 = alns['d1_u2r'].down_coords
    assert (start2, end2) == (9, 14), (start2, end2)
    assert reverse_complement(pair.brk2.upstream[-15:][start2:end2]) == 'GGGTC', pair.brk2.upstream #pair.brk1.upstream[start1:end1]
