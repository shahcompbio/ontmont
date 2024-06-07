import yaml
import pyfaidx
import pandas as pd
import numpy as np

from transposon.irenum import (Breakpoint, BreakpointPair, BreakpointChain, Transitions, Segments)
from transposon.irenum import (get_breakpoint_seqs, enumerate_breakpoints, enumerate_breakpoints, 
    remove_duplicates_from_tra_table, filter_sv_with_breakpoint_at_contig_ends)


CONFIG = """
reference:
    '1': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_GFPVector.fa
    '2': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_PGBD5Vector.fa
    '3': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_PBEF1_PiggyBacVector.fa
    '5': /data1/shahs3/users/chois7/projects/transposon/resources/GRCh38_DelPBEF1_PGBD5Vector.fa
"""
config = yaml.safe_load(CONFIG)


def test_breakpoint_comparison_1():
    brk1 = Breakpoint('chr1', 1000, '+')
    brk2 = Breakpoint('chr1', 2000, '-')
    assert brk1 < brk2, (brk1, brk2)

def test_breakpoint_comparison_2():
    brk1 = Breakpoint('chr1',  1000, '+')
    brk2 = Breakpoint('chr10', 1000, '-')
    assert brk1 < brk2, (brk1, brk2)

def test_breakpoint_comparison_3():
    brk1 = Breakpoint('chr1', 1000, '+')
    brk2 = Breakpoint('chrX', 1000, '-')
    assert brk1 < brk2, (brk1, brk2)

def test_breakpoint_comparison_4():
    brk1 = Breakpoint('chr9',  1000, '+')
    brk2 = Breakpoint('chr10', 1000, '-')
    assert brk1 < brk2, (brk1, brk2)

def test_breakpoint():
    genome = pyfaidx.Fasta(config['reference']['3'])
    brk = Breakpoint('PBEF1NeoTransposon', 1478, '-')
    
    brk.get_breakpoint_seqs(margin=50, genome=genome)
    assert brk.upstream == 'ATGACTGAGCCGGAAAAAAGACCCGACGATATGATCCTGATGCAGCTAGA'
    assert brk.downstream == 'TTAACCCTAGAAAGATAGTCTGCGTAAAATTGACGCATGCATTCTTGAAA'
    assert brk.seq_rearranged == 'TTAACCCTAGAAAGATAGTCTGCGTAAAATTGACGCATGCATTCTTGAAA'
    assert brk.seq_removed == 'ATGACTGAGCCGGAAAAAAGACCCGACGATATGATCCTGATGCAGCTAGA'
    
    brk.get_breakpoint_seqs(margin=10, genome=genome)
    assert brk.upstream == 'TGCAGCTAGA', brk.upstream 
    assert brk.downstream == 'TTAACCCTAG', brk.downstream 
    assert brk.seq_rearranged == 'TTAACCCTAG', brk.seq_rearranged
    assert brk.seq_removed == 'TGCAGCTAGA', brk.seq_removed


def test_get_breakpoint_seqs():
    genome = pyfaidx.Fasta(config['reference']['3'])

    up, dn = get_breakpoint_seqs('chr10', 51340883, margin=50, genome=genome)
    assert up == 'ATAAGACAATCTAATGGGGCAGATGGGGTAAGAGGGTCATTTTATCCATA', up
    assert dn == 'TTAAAAAGATTAAAAACCAAAGTTCTCAGTAATTGATTCTACATCTCACA', dn

    up, dn = get_breakpoint_seqs('PBEF1NeoTransposon', 1478, margin=50, genome=genome)
    assert up == 'ATGACTGAGCCGGAAAAAAGACCCGACGATATGATCCTGATGCAGCTAGA', up
    assert dn == 'TTAACCCTAGAAAGATAGTCTGCGTAAAATTGACGCATGCATTCTTGAAA', dn

    up, dn = get_breakpoint_seqs('PBEF1NeoTransposon', 4995, margin=50, genome=genome)
    assert up == 'CATGATTATCTTTAACGTACGTCACAATATGATTATCTTTCTAGGGTTAA', up
    assert dn == 'TCTAGCTGCGTGTTCTGCAGCGTGTCGAGCATCTTCATCTGCTCCATCAC', dn

    up, dn = get_breakpoint_seqs('PBEF1NeoTransposon', 1, margin=10, genome=genome) # PBEF1NeoTransposon size: 6894
    assert up == 'AAGGATCTTC', up
    assert dn == 'TTGAGATCCT', dn

    up, dn = get_breakpoint_seqs('PBEF1NeoTransposon', 6895, margin=10, genome=genome) # PBEF1NeoTransposon size: 6894
    assert up == 'AAGGATCTTC', up
    assert dn == 'TTGAGATCCT', dn


def test_enumerate_breakpoints():
    df = pd.DataFrame({
        'chrom': {0: 'chr6', 1: 'PBEF1NeoTransposon', 2: 'chr6'},
        'start': {0: 26424060, 1: 4596, 2: 152942012},
        'end': {0: 26424268, 1: 4995, 2: 152944889},
        'strand': {0: '+', 1: '-', 2: '+'},
        'clip1': {0: 30, 1: 2886, 2: 631},
        'match': {0: 209, 1: 398, 2: 2878},
        'clip2': {0: 3280, 1: 235, 2: 10},
        'pclip1': {0: 30, 1: 235, 2: 631}
    })
    brks = enumerate_breakpoints(df)
    assert len(brks) % 2 == 0, brks
    assert brks[0].pos == 26424268, brks[0]
    assert brks[0].ori == '+', brks[0]
    assert brks[1].pos == 4995, brks[1]
    assert brks[1].ori == '+', brks[1]
    assert brks[2].pos == 4596, brks[2]
    assert brks[2].ori == '-', brks[2]
    assert brks[3] == brks[-1], brks[3]
    assert brks[3].pos == 152942012, brks[3]
    assert brks[3].ori == '-', brks[3]

def test_get_transitions():
    df = pd.DataFrame({
        'chrom': {0: 'chr6', 1: 'PBEF1NeoTransposon', 2: 'chr6'},
        'start': {0: 26424060, 1: 4596, 2: 152942012},
        'end': {0: 26424268, 1: 4995, 2: 152944889},
        'strand': {0: '+', 1: '-', 2: '+'},
        'clip1': {0: 30, 1: 2886, 2: 631},
        'match': {0: 209, 1: 398, 2: 2878},
        'clip2': {0: 3280, 1: 235, 2: 10},
        'pclip1': {0: 30, 1: 235, 2: 631}
    })
    brks = enumerate_breakpoints(df)
    brks.get_transitions()
    assert len(brks.tras) == 2, brks.tras
    assert brks.tras[0].brk1.pos == 26424268, brks.tras[0].brk1.pos
    assert brks.tras[0].brk2.pos == 4995, brks.tras[0].brk2.pos
    assert brks.tras[1].brk1.pos == 4596, brks.tras[1].brk1.pos
    assert brks.tras[1].brk2.pos == 152942012, brks.tras[1].brk2.pos

def test_get_segments():
    df = pd.DataFrame({
        'chrom': {0: 'chr6', 1: 'PBEF1NeoTransposon', 2: 'chr6'},
        'start': {0: 26424060, 1: 4596, 2: 152942012},
        'end': {0: 26424268, 1: 4995, 2: 152944889},
        'strand': {0: '+', 1: '-', 2: '+'},
        'clip1': {0: 30, 1: 2886, 2: 631},
        'match': {0: 209, 1: 398, 2: 2878},
        'clip2': {0: 3280, 1: 235, 2: 10},
        'pclip1': {0: 30, 1: 235, 2: 631}
    })
    brks = enumerate_breakpoints(df)
    brks.get_segments()
    assert len(brks.segs) == 1, brks.segs
    assert brks.segs[0].brk1.pos == 4995, brks.segs[0].brk1.pos
    assert brks.segs[0].brk2.pos == 4596, brks.segs[0].brk2.pos


def test_get_transition_list():
    df = pd.DataFrame({
        'qname': {i:'testread' for i in range(4)},
        'chrom': {0: 'chr6', 1: 'PBEF1NeoTransposon', 2: 'chr6', 3: 'chr7'},
        'start': {0: 26424060, 1: 4596, 2: 152942012, 3:52942000},
        'end': {0: 26424268, 1: 4995, 2: 152944889, 3:52950000},
        'strand': {0: '+', 1: '-', 2: '+', 3: '+'},
    })
    transitions = Transitions(df)
    assert len(transitions.list) == 3, transitions.list
    assert transitions.list == [
        (('chr6', 26424268, '+'), ('PBEF1NeoTransposon', 4995, '+')),
        (('PBEF1NeoTransposon', 4596, '-'), ('chr6', 152942012, '-')),
        (('chr6', 152944889, '+'), ('chr7', 52942000, '-')),
    ]


def test_get_segment_list():
    df = pd.DataFrame({
        'qname': {i:'testread' for i in range(4)},
        'chrom': {0: 'chr6', 1: 'PBEF1NeoTransposon', 2: 'chr6', 3: 'chr7'},
        'start': {0: 26424060, 1: 4596, 2: 152942012, 3:52942000},
        'end': {0: 26424268, 1: 4995, 2: 152944889, 3:52950000},
    })
    segments = Segments(df)
    assert segments.df.shape[0] == 4, segments.df
    assert 'qname' in segments.df.columns
    assert segments.list == [
        ('PBEF1NeoTransposon', 4596, 4995), 
        ('chr6', 152942012, 152944889)
    ], segments.list


def test_remove_duplicates_from_tra_tables():
    tra_df_columns = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2', 'support', 
                      'u1_u2_score', 'd1_d2_score', 'u1_d2r_score', 'd1_u2r_score',
                      'u1_u2_pvalue', 'd1_d2_pvalue', 'u1_d2r_pvalue', 'd1_u2r_pvalue']
    data = [
        ['chr1', 100, '+', 'chr2', 200, '-', 1,  0, 0, 3, 4, np.nan, np.nan, 0.01, 0.02],
        ['chr2', 200, '-', 'chr1', 100, '+', 1,  0, 1, 4, 3, np.nan, 0.3, 0.02, 0.01],
    ]
    inputdf = pd.DataFrame(data=data, index=range(0, 2), columns=tra_df_columns)
    outputdf = remove_duplicates_from_tra_table(inputdf)
    
    out_data = [
        ['chr1', 100, '+', 'chr2', 200, '-', 2,  0, 1, 3, 4, np.nan, 0.3, 0.01, 0.02],
    ]
    val = pd.DataFrame(data=out_data, index=range(0, 1), columns=tra_df_columns)
    assert outputdf.shape[0] == val.shape[0], outputdf
    assert np.isnan(outputdf.squeeze()['u1_u2_pvalue']), outputdf.squeeze()['u1_u2_pvalue']
    assert (outputdf.dropna(axis=1) != val.dropna(axis=1)).sum().sum() == 0, outputdf


def test_filter_sv_with_breakpoint_at_contig_ends():
    df = pd.DataFrame({
    0: {
        'chrom': 'PBEF1NeoTransposon',
        'pos1': 1,
        'pos2': 6895,
        'support': 4,
        'segment_score': 9.0,
        'segment_pvalue': 0.0,
    }}).T
    outputdf = filter_sv_with_breakpoint_at_contig_ends(df)
    assert outputdf.shape[0] == 0, outputdf

def test_sort_breakpointchain():
    brks = BreakpointChain([
        Breakpoint('chr2', 1000, '+'),
        Breakpoint('chr1', 1000, '+'),
        Breakpoint('chr1', 2000, '+'),
        Breakpoint('chr1', 1000, '+'),
    ])
    brks.get_transitions(sort_transition=True)
    expected_tra0 = 'chr1:1000:+-chr2:1000:+'
    expected_tra1 = 'chr1:1000:+-chr1:2000:+'
    assert str(brks.tras[0]) == expected_tra0
    assert str(brks.tras[1]) == expected_tra1
