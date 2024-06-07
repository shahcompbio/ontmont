import math
import random
from operator import xor
import yaml
import swalign
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd

import sys
import os

script_dir_path = os.path.abspath(os.path.join(os.path.dirname("__file__"), 'scripts'))
if script_dir_path not in sys.path:
    sys.path.append(script_dir_path)

CONFIG = """
palindrome:
    match: 1
    mismatch: -3
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
"""
config = yaml.load(CONFIG, Loader=yaml.Loader)

scoring_palindrome = swalign.NucleotideScoringMatrix(float(config['palindrome']['match']), float(config['palindrome']['mismatch']))
SW_PALINDROME = swalign.LocalAlignment(scoring_palindrome, globalalign=False,
    gap_penalty=float(config['palindrome']['gap_penalty']), 
    gap_extension_penalty=float(config['palindrome']['gap_extension_penalty']))

scoring_holliday = swalign.NucleotideScoringMatrix(float(config['holliday']['match']), float(config['holliday']['mismatch']))
SW_HOLLIDAY = swalign.LocalAlignment(scoring_holliday, globalalign=False,
    gap_penalty=float(config['holliday']['gap_penalty']), 
    gap_extension_penalty=float(config['holliday']['gap_extension_penalty']))

from label_palindromes import (reverse_complement, 
    label_one_sided_palindromes,  #
    label_two_sided_palindromes,  #
    align_seq_one_sided, align_seq_two_sided, calc_palindrome_prob_one_sided, calc_palindrome_prob_two_sided, align_sequences,
    calc_distance_score, calc_holliday_junction_pval, calc_holliday_junction_stats, is_holliday_junction, label_holliday_junctions, plot_hairpin)

def test_reverse_complement():
    seq_revs = {'ACGGT':'ACCGT', 'GNTTCA': 'TGAA.C'}
    for seq, rev in seq_revs.items():
        assert reverse_complement(seq) == rev, (seq, rev)

def test_align_seq_one_sided():
    _seq1, _seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    seq1, seq2, aln = align_seq_one_sided(_seq1 + _seq2, SW_PALINDROME)
    assert seq1 == _seq1, _seq1
    assert seq2 == _seq2, _seq2
    assert aln.r_pos == 1, aln.dump()
    assert aln.r_end == 13, aln.dump()
    assert aln.q_pos == 1, aln.dump()
    assert aln.q_end == 13, aln.dump()

def test_align_seq_two_sided():
    _seq1, _seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    aln = align_seq_two_sided(_seq1, _seq2, SW_PALINDROME)
    assert aln.r_pos == 1, aln.dump()
    assert aln.r_end == 13, aln.dump()
    assert aln.q_pos == 1, aln.dump()
    assert aln.q_end == 13, aln.dump()

def test_calc_palindrome_prob_one_sided():
    _seq1, _seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    seq1, seq2, aln = align_seq_one_sided(_seq1 + _seq2, SW_PALINDROME)
    assert calc_palindrome_prob_one_sided(seq1 + seq2, SW_PALINDROME, 
        src_score=5, bp_near_breakpoint=4, direction='down', n_iter=100, random_seed=42) == 0.07
    assert calc_palindrome_prob_one_sided(seq1 + seq2, SW_PALINDROME, 
        src_score=6, bp_near_breakpoint=4, direction='down', n_iter=100, random_seed=42) == 0.01
    
def test_calc_palindrome_prob_two_sided():
    _seq1, _seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    seq1, seq2, aln = align_seq_one_sided(_seq1 + _seq2, SW_PALINDROME)
    calc_palindrome_prob_two_sided(seq1, seq2, SW_PALINDROME, 
        src_score=4, n_iter=100, random_seed=42) == 0.36
    calc_palindrome_prob_two_sided(seq1, seq2, SW_PALINDROME, 
        src_score=5, n_iter=100, random_seed=42) == 0.09
    
def test_calc_distance_score():
    assert calc_distance_score(dist1=0, dist2=0, dist_cutoff1=1, dist_cutoff2=3) == 6
    assert calc_distance_score(dist1=1, dist2=3, dist_cutoff1=1, dist_cutoff2=3) == 4
    assert calc_distance_score(dist1=2, dist2=3, dist_cutoff1=1, dist_cutoff2=3) == 2
    assert calc_distance_score(dist1=2, dist2=4, dist_cutoff1=1, dist_cutoff2=3) == 1
    assert calc_distance_score(dist1=4, dist2=4, dist_cutoff1=1, dist_cutoff2=3) == 0

def test_is_holliday_junction():
    score_cutoff = 4
    dist_sum_cutoff = 4 # HARD CODED
    seq1 = 'AACCGGTTACGTACGT'
    seq2 = 'TTTTACGTAACCGGTA'
    negative1, negative2 = False, False
    aln = align_sequences(seq1, seq2, SW_HOLLIDAY)
    assert is_holliday_junction(aln, SW_HOLLIDAY, negative1, negative2, 
        dist_cutoff1=1, dist_cutoff2=9, dist_sum_cutoff=4, score_cutoff=4), aln.dump()

def test_calc_holliday_junction_pval_1():
    # seq1: d1
    # seq2: u2
    # "seq2" (u2 )  TTTTACGT{AACCGGT}A>| (d1 ) | {AACCGGT}TACGTACGT> "seq1"
    #        (u2r) <AAAATGCA{TTGGCCA}T | (d1r) |<{TTGGCCA}ATGCATGCA
    seq1, seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    aln = align_sequences(seq1, seq2, SW_HOLLIDAY)
    assert aln.score == 7, aln.dump()
    assert aln.r_pos == 1-1, aln.dump()
    assert aln.r_end == 7, aln.dump()
    assert aln.orig_ref[1-1:7] == 'AACCGGT', aln.orig_ref
    
    assert aln.q_pos == 9-1, aln.dump()
    assert aln.q_end == 15, aln.dump()
    assert aln.orig_query[9-1:15] == 'AACCGGT', aln.orig_query
    
    sources = ['u5r', 'd5r', 'u5', 'd5']
    xors = [True, False, False, True]
    for src, negative in zip(sources, xors):
        assert xor(src.startswith('d'), src.endswith('r')) == negative
    
    negative1, negative2 = True, False
    dist_cutoff1, dist_cutoff2 = 1, 3
    dist1 = (aln.r_pos) if negative1 else (len(seq1) - aln.r_end)
    dist2 = (aln.q_pos) if negative2 else (len(seq2) - aln.q_end)
    assert dist1 == 0, dist1
    assert dist2 == 1, dist2
    
    dist_score = calc_distance_score(dist1, dist2, dist_cutoff1, dist_cutoff2)
    assert dist_score == (2+1)+(2+1), dist_score
    
    assert calc_holliday_junction_pval(seq1, seq2, negative1, negative2, 
        src_score=5, dist_cutoff1=4, dist_cutoff2=10, sw_holliday=SW_HOLLIDAY,
        n_iter=100, random_seed=42) == 0.08, (seq1, seq2, negative1, negative2)
    
def test_calc_holliday_junction_pval_2():
    # seq1: u1
    # seq2: u2
    # "seq2" (u1 )  {AACCGGT}TACGTACGT>| (u2 )  TTTTACGT{AACCGGT}A>| 
    #        (u1r) <{TTGGCCA}ATGCATGCA | (u2r) <AAAATGCA{TTGGCCA}T | 
    seq1, seq2 = 'AACCGGTTACGTACGT', 'TTTTACGTAACCGGTA'
    aln = align_sequences(seq1, seq2, SW_HOLLIDAY)
    negative1, negative2 = False, False
    dist_cutoff1, dist_cutoff2 = 1, 9
    dist1 = (aln.r_pos) if negative1 else (len(seq1) - aln.r_end)
    dist2 = (aln.q_pos) if negative2 else (len(seq2) - aln.q_end)
    assert dist1 == 9, dist1
    assert dist2 == 1, dist2
    
    dist_score = calc_distance_score(dist1, dist2, dist_cutoff1, dist_cutoff2)
    assert dist_score == (0+1)+(2+1), dist_score
    
    assert calc_holliday_junction_pval(seq1, seq2, negative1, negative2, 
        src_score=5, dist_cutoff1=1, dist_cutoff2=9, sw_holliday=SW_HOLLIDAY,
        n_iter=100, random_seed=42) == 0.03
    
def test_calc_holliday_junction_stats():
    data = pd.DataFrame(data={
        'upstream1': 'AACCGGTTACGTACGT',
        'downstream1': 'CCCCCCCCCCCCCCCC',
        'upstream2': 'TTTTACGTAACCGGTA',
        'downstream2': 'AAAAAAAATTTTTTTT',
    }, index=[0])
    out = data.copy()
    out['holliday_junction_p'] = 0.0
    margin = 16
    row = data.squeeze()
    up1, dn1 = row['upstream1'][-margin:], row['downstream1'][:margin]
    up2, dn2 = row['upstream2'][-margin:], row['downstream2'][:margin]
    min_pval, pvals, alns, seqs = calc_holliday_junction_stats(up1, dn1, up2, dn2, SW_HOLLIDAY, 
        score_cutoff=4, dist_cutoff1=1, dist_cutoff2=9, 
        n_iter=100, random_seed=42, debug=True)
    assert min_pval == 0.0
    assert len(pvals.keys()) == 1 and 'u1_u2' in pvals, pvals
    assert seqs['u1'] == up1, seqs
    assert seqs['u2'] == up2, seqs
    assert seqs['d1'] == dn1, seqs
    assert seqs['d2'] == dn2, seqs
    assert seqs['u2r'] == reverse_complement(up2), seqs
    assert seqs['d2r'] == reverse_complement(dn2), seqs
    assert pvals['u1_u2'] == 0.0, pvals
    assert alns['u1_u2'].score == 7, alns
    assert alns['d1_d2'].score == 0, alns
    assert alns['u1_d2r'].score == 2, alns
    assert alns['d1_u2r'].score == 2, alns

def test_label_holliday_junction():
    data = pd.DataFrame(data={
        'upstream1': 'AACCGGTTACGTACGT',
        'downstream1': 'CCCCCCCCCCCCCCCC',
        'upstream2': 'TTTTACGTAACCGGTA',
        'downstream2': 'AAAAAAAATTTTTTTT',
    }, index=[0])
    out = data.copy()
    out['holliday_junction_p'] = 0.0
    margin = 16
    config = yaml.load(CONFIG, Loader=yaml.Loader)
    config['holliday']['distance_cutoff1'] = 1
    config['holliday']['distance_cutoff2'] = 9
    config['holliday']['score_cutoff'] = 5
    config['holliday']['n_iter'] = 100
    config['holliday']['margin'] = 16
    df = label_holliday_junctions(data, SW_HOLLIDAY, config)
    assert (df != out).sum().sum() == 0, df