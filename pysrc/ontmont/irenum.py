import pandas as pd
import numpy as np
import transposon.irs as irs

class BreakpointChain(list):
    def __init__(self, brks_iterable):
        super().__init__(brks_iterable)
        self.tras = []
        self.segs = []

    # enumeration of tras
    def get_transitions(self, sort_transition=False):
        if len(self) >= 2:
            ix_range = range(0, len(self), 2)
            for i in ix_range:
                brk1 = self[i]
                brk2 = self[i+1]
                if sort_transition:
                    if brk1 > brk2:
                        brk1, brk2 = brk2, brk1
                tra = BreakpointPair(brk1, brk2)
                self.tras.append(tra)

    # enumeration of segs
    def get_segments(self):
        ix_range = range(1, len(self)-1, 2)
        for i in ix_range:
            brk1 = self[i]
            brk2 = self[i+1]
            seg = BreakpointPair(brk1, brk2)
            self.segs.append(seg)


class Sample:
    long_samples = [ # class attribute
        "14472B_201",
        "14472B_202",
        "14472B_500",
        "14472B_501",
        "14472C_203",
        "14472C_204",
        "14472C_502",
        "14472D_101",
        "14472D_102",
        "14472D_103",
        "14472D_104",
        "14472D_205",
        "14472_100",
        "14472_300",
    ]
    sample2group = {s:s.split('_')[1][0] for s in long_samples}
    sample2short = {s:s.split('_')[1] for s in long_samples}
    
    def __init__(self, long_sample):
        self.long_sample = long_sample
        self.group = self.sample2group[self.long_sample]
        self.long_id = self.long_sample
        self.id = self.sample2short[self.long_sample]

    def __repr__(self):
        return f'{self.id}'


class Transitions:
    def __init__(self, df):
        self.df = df.copy()
        self.list = []
        self.get_list()

    def get_list(self):
        for qname, qdf in self.df.groupby('qname'):
            qdf.reset_index(drop=True, inplace=True)
            n_fragments = qdf.shape[0]
            if n_fragments <= 1: continue # don't count TRA when none
            prevdf = qdf.shift(1)
            for rix, row in qdf.iterrows():
                if rix == 0: continue
                prow = prevdf.loc[rix]
                chrom1, start1, end1, strand1 = prow['chrom'], int(prow['start']), int(prow['end']), prow['strand']
                chrom2, start2, end2, strand2 = row['chrom'], int(row['start']), int(row['end']), row['strand']
                if strand1 == '+':
                    ori1 = '+'
                    pos1 = end1
                else:
                    ori1 = '-'
                    pos1 = start1
                if strand2 == '+':
                    ori2 = '-'
                    pos2 = start2
                else:
                    ori2 = '+'
                    pos2 = end2
                tra = ((chrom1, pos1, ori1), (chrom2, pos2, ori2))
                self.list.append(tra)
                

class Segments:
    def __init__(self, df):
        self.df = df.copy()
        self.list = []
        self.get_list()

    def get_list(self):
        for qname, qdf in self.df.groupby('qname'):
            n_fragments = qdf.shape[0]
            if n_fragments <= 2: continue # don't count segments when none
            for rix, row in qdf.iterrows():
                if rix == 0 or rix == n_fragments-1: continue
                chrom, start, end = row['chrom'], int(row['start']), int(row['end'])
                segment = (chrom, start, end)
                self.list.append(segment)


def make_brk_supports(bundle):
    df = pd.DataFrame(columns=['chrom', 'pos', 'ori'])
    for brks in bundle:
        for brk in brks:
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori
            field = [chrom, pos, ori]
            df.loc[df.shape[0]] = field
    brk_supports = df.value_counts()
    return brk_supports


def make_tra_table(bundle, tra_supports):
    holliday_score_cutoff = 5
    key_pairs = ['u1_u2', 'd1_d2', 'u1_d2r', 'd1_u2r']
    tra_cols = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2', 'support']
    tra_cols += [f'{key}_score' for key in key_pairs]
    tra_cols += [f'{key}_pvalue' for key in key_pairs]
    tra_df = pd.DataFrame(columns=tra_cols)
    tra_saved = set()
    for brks in bundle:
        for ix, tra in enumerate(brks.tras):
            coord1 = (tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori)
            coord2 = (tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori)
            coord_pair = (coord1, coord2)
            if coord_pair not in tra_saved:
                tra_saved.add(coord_pair)
            else:
                continue
            support = tra_supports[coord_pair]
            field = [*coord1, *coord2, support]
            field += ([0] * 4)
            field += ([np.nan] * 4)

            for kx, key in enumerate(key_pairs):
                if tra.alns[key]:
                    holliday_score = tra.alns[key].score
                    if holliday_score >= holliday_score_cutoff:
                        field[7+kx] = holliday_score
                        field[11+kx] = tra.alns[key].pvalue
            # if max(field[7:]) > 0:
            tra_df.loc[tra_df.shape[0]] = field
    
    tra_df = filter_sv_with_breakpoint_at_contig_ends(tra_df)
    tra_df = remove_duplicates_from_tra_table(tra_df)
    tra_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return tra_df


def remove_duplicates_from_tra_table(tra_df):
    df = pd.DataFrame(columns = tra_df.columns)
    
    vectors = [
        'DelPBEF1NeoTransposon',
        'GFPVector',
        'PBEF1NeoTransposon',
        'PGBD5Vector',
        'PiggyBacVector',
        'puro-GFP-PGBD5_seq'
    ]
    chrom_order = ['chr'+str(c) for c in range(1, 22)] + ['chrX', 'chrY'] + vectors
    
    for rix, row in tra_df.iterrows():
        chrom1, pos1, ori1 = row['chrom1'], row['pos1'], row['ori1']
        chrom2, pos2, ori2 = row['chrom2'], row['pos2'], row['ori2']
        if chrom1 not in chrom_order or chrom2 not in chrom_order:
            continue
        support = row['support']
        u1_u2_score, d1_d2_score = row['u1_u2_score'], row['d1_d2_score']
        u1_d2r_score, d1_u2r_score = row['u1_d2r_score'], row['d1_u2r_score']
        u1_u2_pvalue, d1_d2_pvalue = row['u1_u2_pvalue'], row['d1_d2_pvalue']
        u1_d2r_pvalue, d1_u2r_pvalue = row['u1_d2r_pvalue'], row['d1_u2r_pvalue']
        flag_reverse = False
        if chrom1 != chrom2:
            if chrom_order.index(chrom1) > chrom_order.index(chrom2):
                flag_reverse = True
        else:
            if pos1 > pos2:
                flag_reverse = True
        if flag_reverse:
            chrom1, pos1, ori1, chrom2, pos2, ori2 = chrom2, pos2, ori2, chrom1, pos1, ori1
            u1_d2r_score, d1_u2r_score = d1_u2r_score, u1_d2r_score
            u1_d2r_pvalue, d1_u2r_pvalue = d1_u2r_pvalue, u1_d2r_pvalue
        field = [chrom1, pos1, ori1, chrom2, pos2, ori2, support, 
            u1_u2_score, d1_d2_score, u1_d2r_score, d1_u2r_score,
            u1_u2_pvalue, d1_d2_pvalue, u1_d2r_pvalue, d1_u2r_pvalue]
        df.loc[df.shape[0]] = field
    
    ix_cols = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2']
    df = df.groupby(ix_cols).agg({
        'support':['sum'], 
        'u1_u2_score':['max'], 'd1_d2_score':['max'], 'u1_d2r_score':['max'], 'd1_u2r_score':['max'],
        'u1_u2_pvalue':['min'], 'd1_d2_pvalue':['min'], 'u1_d2r_pvalue':['min'], 'd1_u2r_pvalue':['min']
    })
    df.columns = df.columns.droplevel(1)
    df = df.reset_index()
    return df


def make_seg_table(bundle, seg_supports):
    segment_score_cutoff = 5
    seg_cols = ['chrom', 'pos1', 'pos2', 'support']
    seg_cols += ['segment_score', 'segment_pvalue']
    seg_df = pd.DataFrame(columns=seg_cols)
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
    
    for brks in bundle:
        for ix, seg in enumerate(brks.segs):
            assert seg.brk1.chrom == seg.brk2.chrom
            chrom = seg.brk1.chrom
            if chrom not in chrom_order: continue #
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
            if seg.aln_segment:
                _segment_score = seg.aln_segment.score
                if _segment_score >= segment_score_cutoff:
                    segment_score = _segment_score
                    segment_pvalue = seg.aln_segment.pvalue
            field = [*coord, support, segment_score, segment_pvalue]
            seg_df.loc[seg_df.shape[0]] = field
    seg_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return seg_df


def make_brk_table(bundle, brk_supports):
    vectors = [
        'DelPBEF1NeoTransposon',
        'GFPVector',
        'PBEF1NeoTransposon',
        'PGBD5Vector',
        'PiggyBacVector',
        'puro-GFP-PGBD5_seq'
    ]
    chrom_order = ['chr'+str(c) for c in range(1, 22)] + ['chrX', 'chrY'] + vectors
    unilateral_score_cutoff = 5
    bilateral_score_cutoff = 8
    brk_cols = ['chrom', 'pos', 'ori', 'support']
    brk_cols += ['upstream_score', 'downstream_score', 'breakpoint_score']
    brk_cols += ['upstream_pvalue', 'downstream_pvalue', 'breakpoint_pvalue']
    brk_df = pd.DataFrame(columns=brk_cols)
    brk_saved = set()
    
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
            field = [*coord, support, 
                upstream_score, downstream_score, breakpoint_score, 
                upstream_pvalue, downstream_pvalue, breakpoint_pvalue]
            brk_df.loc[brk_df.shape[0]] = field

    brk_df = filter_breakpoints_at_contig_ends(brk_df)
    brk_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    
    return brk_df


def make_brks_bundle(df, genome, sw_palindrome, sw_holliday, margins=[15, 30, 60]):    
    bundle = []
    for qname, qdf in df.groupby('qname'):
        brks = enumerate_breakpoints(qdf)
        brks.get_transitions()
        brks.get_segments()
        for brk in brks:
            brk.get_breakpoint_seqs(margin=60, genome=genome)
            seq1 = brk.upstream
            seq2 = brk.downstream
            direction1 = 'up'
            direction2 = 'down'
            brk.aln_upstream = irs.get_best_onesided_ir(seq1, direction1, sw_palindrome, dist_cutoff=2, margins=margins)
            brk.aln_downstream = irs.get_best_onesided_ir(seq2, direction2, sw_palindrome, dist_cutoff=2, margins=margins)
            brk.aln_breakpoint = irs.get_best_ir_within_breakpoints(seq1, seq2, sw_palindrome, dist_cutoff1=100, dist_cutoff2=100, margins=margins)
        for ix, seg in enumerate(brks.segs):
            brks.segs[ix] = irs.get_best_ir_within_segment(
                seg, sw_palindrome, genome, dist_cutoff1=2, dist_cutoff2=5, margins=margins)
        for ix, tra in enumerate(brks.tras):
            brks.tras[ix] = irs.get_best_holliday_junctions(
                tra, sw_holliday, genome, score_cutoff=4, dist_cutoff1=2, dist_cutoff2=5, margins=margins)
        bundle.append(brks)
    return bundle


def enumerate_breakpoints(df):
    df = df.reset_index(drop=True)
    ix_start, ix_end = 0, df.shape[0] - 1
    brks = BreakpointChain([])
    for rix, row in df.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        # match = int(row['match'])
        if rix == ix_start:
            ori = '+' if strand == '+' else '-'
            pos = end if strand == '+' else start
            brk = Breakpoint(chrom, pos, ori)
            _brks = [brk]
        elif rix == ix_end:
            ori = '-' if strand == '+' else '+'
            pos = start if strand == '+' else end
            brk = Breakpoint(chrom, pos, ori)
            _brks = [brk]
        else:
            brk_start = Breakpoint(chrom, start, '-')
            brk_end = Breakpoint(chrom, end, '+')
            _brks = [brk_start, brk_end] if strand == '+' else [brk_end, brk_start]
        brks += _brks
    return brks


def filter_breakpoints_at_contig_ends(brk_df):
    vector_lengths = { # GLOBAL - I don't want to waste time automating this...
        'GFPVector': 8725,
        'PiggyBacVector': 9794,
        'PGBD5Vector': 10218,
        'puro-GFP-PGBD5_seq': 10218,
        'DelPBEF1NeoTransposon': 6796,
        'PBEF1NeoTransposon': 6894,
    }
    cols = brk_df.columns
    brk_df.loc[:, 'remove'] = False
    for rix, row in brk_df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        if chrom not in vector_lengths: continue # only remove circularities
        vector_length = vector_lengths[chrom]
        
        remove = pos <= 5 or abs(pos - vector_length) <= 5
        brk_df.loc[rix, 'remove'] = remove
    brk_df = brk_df[~brk_df['remove']][cols]
    return brk_df


def filter_sv_with_breakpoint_at_contig_ends(df):
    """
    Filters structural variant (SV) records to exclude those with breakpoints too close to the ends of contigs.

    This function calculates the distance of each SV breakpoint from the nearest end of its respective contig.
    SVs with either breakpoint within 3 bases of a contig end are excluded from the results. The function uses
    a predefined dictionary of contig lengths to determine the distances from the contig ends.

    Parameters:
    - df (pandas.DataFrame): A DataFrame containing SV records. Must include 'pos1' and 'chrom1' columns,
      where 'pos1' is the position of the breakpoint and 'chrom1' is the name of the contig.

    Returns:
    - pandas.DataFrame: A filtered DataFrame excluding SVs with breakpoints within 3 bases of the contig ends.
    """
    vector_lengths = { # GLOBAL - I don't want to waste time automating this...
        'GFPVector': 8725,
        'PiggyBacVector': 9794,
        'PGBD5Vector': 10218,
        'puro-GFP-PGBD5_seq': 10218,
        'DelPBEF1NeoTransposon': 6796,
        'PBEF1NeoTransposon': 6894,
    }
    df_cols = df.columns
    for chrom_col in ('chrom1', 'chrom2'):
        if chrom_col not in df_cols:
            assert 'chrom' in df_cols, df_cols
            df[chrom_col] = df['chrom']
            
    data = pd.DataFrame(columns=df_cols)
    df['dist1s'] = np.abs(df['pos1'] - 1)
    df['dist1e'] = np.abs(df['pos1'] - 1 - df['chrom1'].map(vector_lengths))
    df['dist2s'] = np.abs(df['pos2'] - 1)
    df['dist2e'] = np.abs(df['pos2'] - 1 - df['chrom2'].map(vector_lengths))
    df = df.fillna(np.inf)
    df = df[(df['dist1s'] >= 3) & (df['dist1e'] >= 3) & (df['dist2s'] >= 3) & (df['dist2e'] >= 3)]
    df = df[df_cols]
    return df


class Breakpoint:
    chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y', 'M']
    chroms += ['chr'+c for c in chroms]
    def __init__(self, chrom, pos, orientation):
        self.chrom = chrom
        self.pos = pos
        self.ori = orientation
        self.upstream = None
        self.downstream = None
        self.seq_rearranged = None
        self.seq_removed = None

    def get_breakpoint_seqs(self, margin, genome):
        self.upstream, self.downstream = get_breakpoint_seqs(self.chrom, self.pos, margin, genome)
        if self.ori == '+':
            self.seq_rearranged = self.upstream
            self.seq_removed = self.downstream
        elif self.ori == '-':
            self.seq_rearranged = self.downstream
            self.seq_removed = self.upstream
        else:
            raise ValueError(f'self.ori = {self.ori}')

    def __repr__(self):
        return f'{self.chrom}:{self.pos}:{self.ori}'

    def __lt__(self, other):
        self_chrom_ix = self.chroms.index(self.chrom)
        other_chrom_ix = self.chroms.index(other.chrom)
        if self_chrom_ix < other_chrom_ix:
            return True
        elif self_chrom_ix == other_chrom_ix and self.pos < other.pos:
            return True
        return False


def get_breakpoint_seqs(chrom, pos, margin, genome):
    linear_chroms = set(['chr'+str(c) for c in range(1, 22+1)] + ['chrX', 'chrY'])
    vector_lengths = { # GLOBAL - I don't want to waste time automating this...
        'GFPVector': 8725,
        'PiggyBacVector': 9794,
        'PGBD5Vector': 10218,
        'puro-GFP-PGBD5_seq': 10218,
        'DelPBEF1NeoTransposon': 6796,
        'PBEF1NeoTransposon': 6894,
    }
    assert chrom in genome, chrom
    # assert pos in genome[chrom], pos
    zpos = pos-1 # zero-based
    chrom_len = len(genome[chrom])
    start = zpos - margin
    end = zpos + margin
    if chrom in linear_chroms:
        start = max(0, start)
        end = min(chrom_len, end)
    
    if start < 0:
        upstream = genome[chrom][start:].seq.upper() + genome[chrom][:zpos].seq.upper() # upstream of 0
    else:
        upstream = genome[chrom][start : zpos].seq.upper() # upstream
    if end >= chrom_len:
        residual = end - chrom_len
        downstream = genome[chrom][zpos:].seq.upper() + genome[chrom][:residual].seq.upper() # downstream over end
    else:
        downstream = genome[chrom][zpos : end].seq.upper() # downstream
    return upstream, downstream

class BreakpointPair:
    def __init__(self, brk1, brk2):
        self.brk1 = brk1
        self.brk2 = brk2

    def __repr__(self):
        return f'{self.brk1.chrom}:{self.brk1.pos}:{self.brk1.ori}-{self.brk2.chrom}:{self.brk2.pos}:{self.brk2.ori}'

