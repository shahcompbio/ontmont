import re
import pandas as pd
import numpy as np


class SplitAlignment:
    """
    Represents a split alignment of a sequencing read, encapsulating details about the alignment's location,
    orientation, and structure as described by its CIGAR string.

    This class is designed to parse and store information about a single alignment segment of a read, which
    can be part of a primary alignment or a secondary (chimeric) alignment. It includes methods to interpret
    the CIGAR string for further analysis of the alignment's properties, such as calculating the aligned
    length and identifying specific alignment features.
    """
    def __init__(self, cigarstring, read_name, refname, read_pos, strand):
        """
        Attributes:
        - cigarstring (str): The CIGAR string representing the alignment's match, insertion, deletion, 
            skipping, soft clipping, hard clipping, padding, sequence match, and sequence mismatch operations.
        - read_name (str): The name (query name) of the read.
        - refname (str): The name of the reference chromosome to which the read is aligned.
        - read_pos (int): The starting position of the alignment on the reference chromosome (0-based).
        - strand (str): The strand of the alignment, represented as '+' for the forward strand and '-' for
            the reverse strand.
        """
        self.read_name = read_name
        self.refname = refname
        self.cigarstring = cigarstring
        self.cigar_tuples = SplitAlignment.get_cigar_tuples(self.cigarstring)
        self.primary = None
        self.start = read_pos
        self.strand = strand
        self.match = 0

        self.aln_cols = ['qname', 'chrom', 'start', 'end', 'strand', 'clip1', 'match', 'clip2', 'pclip1']
        self.extract_cigar_field()
        
    def extract_cigar_field(self):
        split_ctags = {4, 5}
        query_consumers = {0, 1, 4, 7, 8}
        reference_consumers = {0, 2, 3, 7, 8}
        self.clip1, self.match, self.clip2 = 0, 0, 0
        flag_after_match = False
        chrom, pos, strand, cigar_tuples = self.refname, self.start, self.strand, self.cigar_tuples
        tpos = int(pos) # init reference pos
        for ctag, clen in cigar_tuples:
            if ctag in split_ctags:
                if not flag_after_match:
                    self.clip1 += clen
                else:
                    self.clip2 += clen
            elif ctag in query_consumers:
                flag_after_match = True
                self.match += clen
            if ctag in reference_consumers:
                flag_after_match = True
                tpos += clen
                continue
        self.end = tpos
        pclip1 = self.clip1 if strand == '+' else self.clip2
        self.pclip1 = pclip1 # positively corrected first clip length

    @staticmethod
    def get_cigar_tuples(cigarstring):
        """
        Returns a cigar tuple from a CIGAR string
        """
        cigar_types = 'MIDNSHP=X'
        converter = {
            cigar_types[i]:i for i in range(len(cigar_types))
        }
        cigar_tuples = []
        for match in re.finditer('(\d+)([A-Z\=])', cigarstring):
            clen, ctype = match.groups()
            cigar_tuple = (converter[ctype], int(clen))
            cigar_tuples.append(cigar_tuple)
        return cigar_tuples


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


class BreakpointPair:
    def __init__(self, brk1, brk2):
        self.brk1 = brk1
        self.brk2 = brk2

    def __repr__(self):
        return f'{self.brk1.chrom}:{self.brk1.pos}:{self.brk1.ori}-{self.brk2.chrom}:{self.brk2.pos}:{self.brk2.ori}'

