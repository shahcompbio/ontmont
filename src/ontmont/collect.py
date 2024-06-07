from collections import Counter

import tqdm
import pandas as pd
import pysam

from .utils import get_secondaries, make_split_read_table, enumerate_breakpoints
from .datatypes import SplitAlignment, Breakpoint, BreakpointPair

def extract_split_alignments(reads, max_depth=500):
    alignments = []
    for i, read in enumerate(reads):
        if i >= max_depth:
            break
        secondary_list = get_secondaries(read)
        if len(secondary_list) < 1: # only consider reads with secondaries (chimeric reads)
            continue # return None
        strand = '-' if read.is_reverse else '+'
        tpos = read.pos + 1
        cigarstring = read.cigarstring
        qname = read.qname
        refname = read.reference_name
        alignment = SplitAlignment(cigarstring, qname, refname, tpos, strand)
        alignments.append(alignment)

        for _, secondary in enumerate(secondary_list):
            refname_s, tpos_s, strand_s, cigarstring_s, _, _ = secondary.split(',')
            tpos_s = int(tpos_s)
            alignment = SplitAlignment(cigarstring_s, qname, refname_s, tpos_s, strand_s)
            alignments.append(alignment)
    return alignments

def pull_breakpoints_from_reads_in_sv_regions(bam, tra, get_read_table=False, min_n_breakpoint=3, margin=10):
    complexes = []
    saved_qnames = set()
    canonical_chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y']
    canonical_chroms += ['chr'+c for c in canonical_chroms]
    read_info = pd.DataFrame()
    for tix, row in tqdm.tqdm(tra.iterrows(), total=tra.shape[0]):
        # row = tra.iloc[6] # breakpoint pair: SV
        chrom1 = row['chromosome_1']
        pos1 = row['position_1']
        chrom2 = row['chromosome_2']
        pos2 = row['position_2']
        reads1 = bam.fetch(chrom1, pos1-margin-1, pos1+margin)
        reads2 = bam.fetch(chrom2, pos2-margin-1, pos2+margin)
        for reads in [reads1, reads2]:
            alignments = extract_split_alignments(reads)
            df = make_split_read_table(alignments)
            df = df[~df['qname'].isin(saved_qnames)] # don't double-count breakpoints already saved
            read_info = pd.concat([read_info, df])
            df = df[df['chrom'].isin(canonical_chroms)] # remove segments in decoys
            _qnames = set(df['qname'].unique().tolist())
            saved_qnames.update(_qnames)
            bundle = make_brks_bundle(df)
            for brks in bundle:
                if len(brks) >= min_n_breakpoint:
                    complexes.append(brks)
    if get_read_table:
        return complexes, read_info
    return complexes

def make_brks_bundle(df):    
    bundle = []
    for qname, qdf in df.groupby('qname'):
        brks = enumerate_breakpoints(qdf)
        brks.qname = qname
        brks.info = {'sv':brks.tras, 'qname':brks.qname}
        brks.get_transitions(sort_transition=False)
        bundle.append(brks)
    return bundle

def get_breakpoint_support_from_bundle(complexes):
    breakpoint_support = Counter()
    for brks in complexes:
        for brk in brks:
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori
            breakpoint_support[str(brk)] += 1
    return breakpoint_support

def map_similar_coordinate_to_higher_rank(complexes, breakpoint_support, margin=10):
    coord_map = {}
    coord_map_log = {}
    for cix, brks in enumerate(complexes):
        for bix, brk in enumerate(brks):
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori 
            src_coord = str(brk)
            src_count = breakpoint_support[src_coord]
            max_count = src_count
            max_coord = src_coord
            for _pos in range(pos-margin, pos+margin):
                _coord = f'{chrom}:{_pos}:{ori}'
                if _coord in breakpoint_support:
                    _count = breakpoint_support[_coord]
                    if _count > max_count:
                        max_count = _count
                        max_coord = _coord
            coord_map[src_coord] = max_coord
            coord_map_log[src_coord] = (max_coord, src_count, max_count)
    return coord_map, coord_map_log

def fix_lower_support_coordinates(complexes, coord_map):
    for cix, brks in enumerate(complexes):
        for bix, brk in enumerate(brks):
            src_coord = str(brk)
            dst_coord = coord_map[src_coord]
            dst_chrom, dst_pos, dst_ori = dst_coord.split(':')
            dst_pos = int(dst_pos)
            complexes[cix][bix].pos = dst_pos
    return complexes

def normalize_sv_table(df, chrom1_col='chromosome_1', chrom2_col='chromosome_2', 
                      pos1_col='position_1', pos2_col='position_2', 
                      ori1_col='strand_1', ori2_col='strand_2'):
    df = df.copy()
    chroms = [str(c) for c in range(1, 22+1)]+['X', 'Y']
    chroms += ['chr'+c for c in chroms] # chr prefices
    chrom_map = dict(zip(chroms, range(len(chroms))))
    assert (~df[chrom1_col].isin(chroms)).sum() == 0, df[chrom1_col].unique()
    assert (~df[chrom2_col].isin(chroms)).sum() == 0, df[chrom2_col].unique()
    diff_chr = df[chrom1_col] != df[chrom2_col]
    flag_inverse = (df[chrom1_col].map(chrom_map) > df[chrom2_col].map(chrom_map))
    flag_inverse |= (df[chrom1_col]==df[chrom2_col]) & (df[pos1_col] > df[pos2_col])
    df.loc[flag_inverse, [chrom1_col, chrom2_col]] = df.loc[flag_inverse, [chrom2_col, chrom1_col]].values
    df.loc[flag_inverse, [pos1_col, pos2_col]] = df.loc[flag_inverse, [pos2_col, pos1_col]].values
    df.loc[flag_inverse, [ori1_col, ori2_col]] = df.loc[flag_inverse, [ori2_col, ori1_col]].values
    return df

def get_svtype(tra:BreakpointPair):
    translocation = 'TRA'
    inversion = 'INV'
    duplication = 'DUP'
    insertion = 'INS'
    deletion = 'DEL'
    chrom1 = tra.brk1.chrom
    chrom2 = tra.brk2.chrom
    if chrom1 != chrom2:
        return translocation
    else: # same chrom
        pos1, pos2 = tra.brk1.pos, tra.brk2.pos
        ori1, ori2 = tra.brk1.ori, tra.brk2.ori
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
            ori1, ori2 = ori2, ori1
        if ori1 == ori2:
            return inversion
        else:
            if ori1 == '+':
                if ori2 == '-':
                    return deletion
            elif ori1 == '-':
                if ori2 == '+':
                    return duplication
    raise ValueError(f'tra:{tra}')

def pull_sv_supporting_reads_from_bundle(sv, bundle):
    output_bundle = []
    sv_prop = {
        (sv['chromosome_1'], sv['position_1'], sv['strand_1']), 
        (sv['chromosome_2'], sv['position_2'], sv['strand_2']), 
    }
    for brks in bundle:
        brks.contains_sv_support = False
        for tra in brks.tras:
            tra.supports_sv = False
            brk1 = tra.brk1
            brk2 = tra.brk2
            tra_prop = {
                (brk1.chrom, brk1.pos, brk1.ori),
                (brk2.chrom, brk2.pos, brk2.ori),
            }
            if sv_prop == tra_prop:
                tra.supports_sv = True
                brks.contains_sv_support = True
        if brks.contains_sv_support:
            output_bundle.append(brks)
    return output_bundle

def find_presence_of_matching_sv(sv1, sv2, margin=50):
    match = []
    for i, (rix, row) in enumerate(sv2.iterrows()):
        _match = (
            (sv1['chromosome_1'] == row['chromosome_1']) &
            (sv1['chromosome_2'] == row['chromosome_2']) &
            (abs(sv1['position_1']-row['position_1']) < margin) &
            (abs(sv1['position_2']-row['position_2']) < margin) &
            (sv1['strand_1'] == row['strand_1']) &
            (sv1['strand_2'] == row['strand_2'])
        )
        if i == 0:
            match = _match
        else:
            match |= _match
    return match

def pull_breakpoints_from_bam_files(bam_paths, sv, get_read_table=False):
    complexes = []
    read_tables = pd.DataFrame()
    for bam_path in bam_paths:
        bam = pysam.AlignmentFile(bam_path)
        if get_read_table:
            _complex, read_table = pull_breakpoints_from_reads_in_sv_regions(bam, sv, get_read_table=get_read_table)
            read_tables = pd.concat([read_tables, read_table])
        else:
            _complex = pull_breakpoints_from_reads_in_sv_regions(bam, sv, get_read_table=get_read_table)
        complexes += _complex
    if get_read_table:
        return complexes, read_tables
    return complexes

def make_tumor_sv_table(complexes, sv, margin=10, get_support=True):
    breakpoint_support = get_breakpoint_support_from_bundle(complexes)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(complexes, breakpoint_support, margin=margin)
    complexes = fix_lower_support_coordinates(complexes, coord_map)
    # breakpoint_support_post = get_breakpoint_support_from_bundle(complexes)
    sv_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    tumor_sv = pd.DataFrame(columns=sv_cols)
    for brks in complexes:
        for tra in brks.tras:
            svtype = get_svtype(tra)
            chrom1, chrom2 = tra.brk1.chrom, tra.brk2.chrom
            pos1, pos2 = tra.brk1.pos, tra.brk2.pos
            ori1, ori2 = tra.brk1.ori, tra.brk2.ori
            field = [chrom1, pos1, ori1, chrom2, pos2, ori2, svtype]
            tumor_sv.loc[tumor_sv.shape[0]] = field
    tumor_sv = normalize_sv_table(tumor_sv)
    ix_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    if get_support:
        tumor_sv = tumor_sv.groupby(ix_cols).value_counts().reset_index()
        tumor_sv.rename(columns={'count':'support'}, inplace=True)
    tumor_sv['in_source'] = find_presence_of_matching_sv(tumor_sv, sv, margin=50)
    return tumor_sv

def get_normalized_sv(tra):
    chrom1, pos1, ori1 = tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori
    chrom2, pos2, ori2 = tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori
    if Breakpoint.chroms.index(chrom1) < Breakpoint.chroms.index(chrom1):
        chrom1, pos1, ori1, chrom2, pos2, ori2 = chrom2, pos2, ori2, chrom1, pos1, ori1
    elif chrom1 == chrom2:
        if pos1 > pos2:
            chrom1, pos1, ori1, chrom2, pos2, ori2 = chrom2, pos2, ori2, chrom1, pos1, ori1
    return [chrom1, pos1, ori1, chrom2, pos2, ori2]

