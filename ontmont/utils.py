import random

import pandas as pd
import numpy as np

from .datatypes import Breakpoint, BreakpointChain, SplitAlignment

def remove_duplicates_from_tra_table(tra_df):
    data = []
    
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
        data.append(field)
    df = pd.DataFrame(data, columns=tra_df.columns)
    
    ix_cols = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2']
    df = df.groupby(ix_cols).agg({
        'support':['sum'], 
        'u1_u2_score':['max'], 'd1_d2_score':['max'], 'u1_d2r_score':['max'], 'd1_u2r_score':['max'],
        'u1_u2_pvalue':['min'], 'd1_d2_pvalue':['min'], 'u1_d2r_pvalue':['min'], 'd1_u2r_pvalue':['min']
    })
    df.columns = df.columns.droplevel(1)
    df = df.reset_index()
    return df


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
            
    df['dist1s'] = np.abs(df['pos1'] - 1)
    df['dist1e'] = np.abs(df['pos1'] - 1 - df['chrom1'].map(vector_lengths))
    df['dist2s'] = np.abs(df['pos2'] - 1)
    df['dist2e'] = np.abs(df['pos2'] - 1 - df['chrom2'].map(vector_lengths))
    df[['dist1s', 'dist1e', 'dist2s', 'dist2e']] = df[['dist1s', 'dist1e', 'dist2s', 'dist2e']].astype(float).fillna(np.inf)
    df = df[(df['dist1s'] >= 3) & (df['dist1e'] >= 3) & (df['dist2s'] >= 3) & (df['dist2e'] >= 3)]
    df = df[df_cols]
    return df


def get_secondaries(read):
    secondaries = [a[1].split(';') for a in read.tags if a[0] == 'SA']
    if len(secondaries) > 0:
        secondaries = secondaries[0]
        secondaries = list(filter(lambda x: x!='', secondaries))
        return secondaries
    return []


def get_chromosomes_to_process(bam, drop_expression_vectors=False):
    vector_contigs = ['GFPVector', 'PiggyBacVector', 'PGBD5Vector', 
        'puro-GFP-PGBD5_seq', 'DelPBEF1NeoTransposon', 'PBEF1NeoTransposon']
    if drop_expression_vectors:
        vector_contigs = ['DelPBEF1NeoTransposon', 'PBEF1NeoTransposon']
    chroms_proc = []
    for contig in bam.references:
        if contig not in vector_contigs: continue
        if bam.count(contig=contig):
            print(contig)
            chroms_proc.append(contig)
    return chroms_proc


def extract_split_alignments(bam, chroms_proc):
    alignments = []
    for chrom in chroms_proc:
        for read in bam.fetch(chrom):
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

def make_split_read_table(alignments):
    sa_cols = ['qname', 'chrom', 'start', 'end', 'strand', 'clip1', 'match', 'clip2', 'pclip1']
    data = []
    for alignment in alignments:
        qname = alignment.read_name
        chrom = alignment.refname
        start = alignment.start
        end = alignment.end
        strand = alignment.strand
        clip1 = alignment.clip1
        clip2 = alignment.clip2
        match = alignment.match
        pclip1 = alignment.pclip1
        field = [qname, chrom, start, end, strand, clip1, match, clip2, pclip1]
        data.append(field)
        # df.loc[df.shape[0]] = field
    df = pd.DataFrame(data, columns=sa_cols)
    df.drop_duplicates(inplace=True)
    df.sort_values(by=['qname', 'pclip1'], inplace=True)
    return df.reset_index(drop=True)

def is_breakpoints_not_sorted(chrom1, pos1, chrom2, pos2, chrom_order):
    if chrom_order.index(chrom1) > chrom_order.index(chrom2):
        return True
    if chrom_order.index(chrom1) == chrom_order.index(chrom2) and pos1 > pos2:
        return True
    return False

def reverse_complement(seq):
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'.'}
    revcmp = ''.join([comp[b] for b in seq[::-1]])
    return revcmp

def shuffle_seq(seq):
    seq_list = [s for s in seq]
    random.shuffle(seq_list)
    shuffled_seq = ''.join(seq_list)
    return shuffled_seq
