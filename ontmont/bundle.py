import pandas as pd
import numpy as np
import tqdm

from .utils import (remove_duplicates_from_tra_table, 
    filter_sv_with_breakpoint_at_contig_ends, filter_breakpoints_at_contig_ends, enumerate_breakpoints)
from .irs import (get_best_onesided_ir, get_best_ir_within_breakpoints, 
    get_best_holliday_junctions, get_best_ir_within_segment)

def make_seg_table(bundle, seg_supports, segment_score_cutoff=5, label_irs=False):
    """Create a dataframe based on a ``BreakpointChain`` bundle and supports dict

    Args:
        bundle (list): List of ``BreakpointChain`` variables
        seg_supports (dict | pandas.Series): Number of support for each breakpoint coordinate
        segment_score_cutoff (int, optional): Alignment score cutoff the IR found in the segment. Defaults to 5.

    Returns:
        pandas.DataFrame: table of segments coordinate with supports and IR statistics
    """
    seg_saved = set()
    vectors = [
        'DelPBEF1NeoTransposon',
        'GFPVector',
        'PBEF1NeoTransposon',
        'PGBD5Vector',
        'PiggyBacVector',
        'puro-GFP-PGBD5_seq'
    ]
    chrom_order = ['chr'+str(c) for c in range(1, 22)] + ['chrX', 'chrY'] + vectors
    data = []
    
    for brks in bundle:
        for seg in brks.segs:
            assert seg.brk1.chrom == seg.brk2.chrom
            chrom = seg.brk1.chrom
            if chrom not in chrom_order: 
                continue
            pos1 = seg.brk1.pos
            pos2 = seg.brk2.pos
            ori1 = seg.brk1.ori
            ori2 = seg.brk2.ori
            if pos1 > pos2:
                pos1, ori1, pos2, ori2 = pos2, ori2, pos1, ori1
            coord = (chrom, pos1, pos2)
            if coord not in seg_saved:
                seg_saved.add(coord)
            else:
                continue
            segment_score = 0
            segment_pvalue = np.nan
            if coord not in seg_supports:
                continue
            support = seg_supports[coord]
            field = [*coord, support]

            if seg.aln_segment:
                _segment_score = seg.aln_segment.score
                if _segment_score >= segment_score_cutoff:
                    segment_score = _segment_score
                    segment_pvalue = seg.aln_segment.pvalue
            
            if label_irs:
                field += [segment_score, segment_pvalue]
            data.append(field)
    seg_cols = ['chrom', 'pos1', 'pos2', 'support']
    if label_irs:
        seg_cols += ['segment_score', 'segment_pvalue']
    seg_df = pd.DataFrame(data, columns=seg_cols)
    seg_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return seg_df


def make_brk_table(bundle, brk_supports,
    unilateral_score_cutoff=5, bilateral_score_cutoff=8, label_irs=False):
    """Create a dataframe of breakpoints

    Args:
        bundle (list): List of ``BreakpointChain`` variables
        brk_supports (dict | pandas.Series): Number of support for each breakpoint coordinate
        unilateral_score_cutoff (int, optional): Alignment score cutoff for unilateral IR. Defaults to 5.
        bilateral_score_cutoff (int, optional): Alignment score cutoff for bilateral IR. Defaults to 8.

    Returns:
        pandas.DataFrame: table of breakpoint coordinates with support and IR statistics
    """
    vectors = [
        'DelPBEF1NeoTransposon',
        'GFPVector',
        'PBEF1NeoTransposon',
        'PGBD5Vector',
        'PiggyBacVector',
        'puro-GFP-PGBD5_seq'
    ]
    chrom_order = ['chr'+str(c) for c in range(1, 22)] + ['chrX', 'chrY'] + vectors
    brk_saved = set()
    data = []
    for brks in bundle:
        for brk in brks:
            if brk.chrom not in chrom_order: continue
            coord = (brk.chrom, brk.pos, brk.ori)
            if coord not in brk_saved:
                brk_saved.add(coord)
            else:
                continue
            support = brk_supports[coord]
            upstream_score, downstream_score, breakpoint_score = 0, 0, 0
            upstream_pvalue, downstream_pvalue, breakpoint_pvalue = np.nan, np.nan, np.nan
            
            if hasattr(brk, 'aln_upstream') and brk.aln_upstream:
                _upstream_score = brk.aln_upstream.score
                if _upstream_score > unilateral_score_cutoff:
                    upstream_score = _upstream_score
                    upstream_pvalue = brk.aln_upstream.pvalue
                    
            if hasattr(brk, 'aln_downatream') and brk.aln_downatream:
                _downstream_score = brk.aln_downstream.score
                if _downstream_score > unilateral_score_cutoff:
                    downstream_score = _downstream_score
                    downstream_pvalue = brk.aln_downstream.pvalue
                    
            if hasattr(brk, 'aln_breakpoint') and brk.aln_breakpoint:
                _breakpoint_score = brk.aln_breakpoint.score
                if _breakpoint_score >= bilateral_score_cutoff:
                    breakpoint_score = _breakpoint_score
                    breakpoint_pvalue = brk.aln_breakpoint.pvalue
    
            # if max([upstream_score, downstream_score, breakpoint_score]) > 0:
            field = [*coord, support]
            if label_irs:
                field += [upstream_score, downstream_score, breakpoint_score, 
                          upstream_pvalue, downstream_pvalue, breakpoint_pvalue]
            data.append(field)

    brk_cols = ['chrom', 'pos', 'ori', 'support']
    if label_irs:
        brk_cols += ['upstream_score', 'downstream_score', 'breakpoint_score']
        brk_cols += ['upstream_pvalue', 'downstream_pvalue', 'breakpoint_pvalue']
    brk_df = pd.DataFrame(data, columns=brk_cols)
    brk_df = filter_breakpoints_at_contig_ends(brk_df)
    brk_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    
    return brk_df


def make_aligned_brks_bundle(reads_df, genome=None, sw_palindrome=None, sw_holliday=None, margins=(15, 30, 60), track_progress=False):
    """Make a list of ``BreapointChain`` based on alignment table, genome, and alignment parameters

    Args:
        reads_df (pandas.DataFrame): Table of read alignment statistics
        genome (pyfaidx.Fasta): Genome fasta
        sw_palindrome (swalign.LocalAlignment): Parameters for detecting IR
        sw_holliday (swalign.LocalAlignment): Parameters for detecting homology
        margins (list, optional): Bases to slice from breakpoints. Defaults to (15, 30, 60).

    Returns:
        list: List of ``BreakpointChain``
    """
    bundle = []
    margin_max = max(margins)
    iterable = reads_df.groupby('qname')
    iterable = tqdm.tqdm(iterable) if track_progress else iterable
    for qname, qdf in iterable:
        brks = enumerate_breakpoints(qdf)
        brks.qname = qname
        brks.get_transitions()
        brks.get_segments()

        flag_align = sw_palindrome and sw_holliday
        if flag_align: # find IRs through SW alignment
            for brk in brks:
                brk.get_breakpoint_seqs(margin=margin_max, genome=genome)
                seq1 = brk.upstream
                seq2 = brk.downstream
                direction1 = 'up'
                direction2 = 'down'
                brk.aln_upstream = get_best_onesided_ir(seq1, direction1, sw_palindrome, dist_cutoff=2, margins=margins)
                brk.aln_downstream = get_best_onesided_ir(seq2, direction2, sw_palindrome, dist_cutoff=2, margins=margins)
                brk.aln_breakpoint = get_best_ir_within_breakpoints(seq1, seq2, sw_palindrome, dist_cutoff1=100, dist_cutoff2=100, margins=margins)
            for ix, seg in enumerate(brks.segs):
                brks.segs[ix] = get_best_ir_within_segment(
                    seg, sw_palindrome, genome, dist_cutoff1=2, dist_cutoff2=5, margins=margins)
            for ix, tra in enumerate(brks.tras):
                brks.tras[ix] = get_best_holliday_junctions(
                    tra, sw_holliday, genome, score_cutoff=4, dist_cutoff1=2, dist_cutoff2=5, margins=margins)

        bundle.append(brks)
    return bundle


def make_brk_supports(bundle):
    """Count supports for unique breakpoint coordinates

    Args:
        bundle (list): List of ``BreakpointChain``

    Returns:
        pandas.Series: Count of unique breakpoints, taking chrom, pos, ori into account
    """
    df = pd.DataFrame(columns=['chrom', 'pos', 'ori'])
    for brks in bundle:
        for brk in brks:
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori
            field = [chrom, pos, ori]
            df.loc[df.shape[0]] = field
    brk_supports = df.value_counts()
    return brk_supports


def make_tra_table(bundle, tra_supports, label_irs=False):
    """Make a table of SVs based on bundle and number of supports

    Args:
        bundle (list): List of ``BreakpointChain``
        tra_supports (dict | pandas.Series): Support count for SVs

    Returns:
        pandas.DataFrame: Table of SVs, duplicate removed, removed if located on contig termini
    """
    holliday_score_cutoff = 5
    key_pairs = ['u1_u2', 'd1_d2', 'u1_d2r', 'd1_u2r']
    tra_saved = set()
    data = []
    for brks in bundle:
        for tra in brks.tras:
            coord1 = (tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori)
            coord2 = (tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori)
            coord_pair = (coord1, coord2)
            if coord_pair not in tra_saved:
                tra_saved.add(coord_pair)
            else:
                continue
            support = tra_supports[coord_pair]
            field = [*coord1, *coord2, support]
            if label_irs:
                field += ([0] * 4) # add aln scores
                field += ([np.nan] * 4) # add p values

                for kx, key in enumerate(key_pairs):
                    if tra.alns[key]:
                        holliday_score = tra.alns[key].score
                        if holliday_score >= holliday_score_cutoff:
                            field[7+kx] = holliday_score
                            field[11+kx] = tra.alns[key].pvalue
            # if max(field[7:]) > 0:
            data.append(field)
    
    tra_cols = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2', 'support']
    if label_irs:
        tra_cols += [f'{key}_score' for key in key_pairs]
        tra_cols += [f'{key}_pvalue' for key in key_pairs]
    tra_df = pd.DataFrame(data, columns=tra_cols)
    tra_df = filter_sv_with_breakpoint_at_contig_ends(tra_df)
    if label_irs:
        tra_df = remove_duplicates_from_tra_table(tra_df)
    tra_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return tra_df
