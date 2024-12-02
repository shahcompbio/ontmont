import os
import pytest

import pysam
import numpy as np
import pandas as pd

from ontmont.datatypes import Breakpoint, BreakpointPair, BreakpointChain

from ontmont.collect import (map_similar_coordinate_to_higher_rank, fix_lower_support_coordinates, 
    get_breakpoint_support_from_bundle, normalize_sv_table, pull_sv_supporting_reads_from_bundle, find_presence_of_matching_sv,
    extract_read_data)

@pytest.fixture
def bundle():
    brk_chain1 = [
        Breakpoint('chr1', 1000, '+'),
        Breakpoint('chr2', 2000, '-'),
        Breakpoint('chr1', 1000, '-'),
        Breakpoint('chr3', 1000, '+'),
    ]
    brk_chain2 = [
        Breakpoint('chr1', 1001, '+'),
        Breakpoint('chr1', 1011, '+'),
    ]
    brk_chain3 = [
        Breakpoint('chr1', 1001, '+'),
        Breakpoint('chr1', 1012, '+'),
    ]
    bundle = [BreakpointChain(brk_chain1), BreakpointChain(brk_chain2), BreakpointChain(brk_chain3)]
    for brk_chain in bundle:
        brk_chain.get_transitions()
    return bundle

def test_map_similar_coordinate_to_higher_rank(bundle):
    #bundle # fixture
    breakpoint_support = get_breakpoint_support_from_bundle(bundle)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=10)
    assert coord_map['chr1:1000:+'] == 'chr1:1001:+', coord_map
    assert coord_map['chr1:1000:-'] == 'chr1:1000:-', coord_map
    assert coord_map['chr1:1011:+'] == 'chr1:1001:+', coord_map
    assert coord_map['chr1:1012:+'] == 'chr1:1012:+', coord_map

def test_fix_lower_support_coordinates(bundle):
    breakpoint_support = get_breakpoint_support_from_bundle(bundle)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=10)
    fixed_bundle = fix_lower_support_coordinates(bundle, coord_map)
    expected_output1 = BreakpointPair(Breakpoint('chr1', 1001, '+'), Breakpoint('chr2', 2000, '-')) # 1000 fixed to 1001
    expected_output2 = BreakpointPair(Breakpoint('chr1', 1001, '+'), Breakpoint('chr1', 1001, '+')) # 1011 fixed to 1001
    assert str(fixed_bundle[0].tras[0]) == str(expected_output1), fixed_bundle[0].tras[0]
    assert str(fixed_bundle[1].tras[0]) == str(expected_output2), fixed_bundle[1].tras[0]

def test_normalize_sv_table():
    input_sv = pd.DataFrame({
    0: {
        'chromosome_1': '11',
        'position_1': 2000,
        'strand_1': '+',
        'chromosome_2': '2',
        'position_2': 1000,
        'strand_2': '-'
    },
    1: {
        'chromosome_1': '1',
        'position_1': 2000,
        'strand_1': '+',
        'chromosome_2': '1',
        'position_2': 1000,
        'strand_2': '-'
    }}).T
    expected_sv = pd.DataFrame({
    0: {
        'chromosome_1': '2',
        'position_1': 1000,
        'strand_1': '-',
        'chromosome_2': '11',
        'position_2': 2000,
        'strand_2': '+'
    },
    1: {
        'chromosome_1': '1',
        'position_1': 1000,
        'strand_1': '-',
        'chromosome_2': '1',
        'position_2': 2000,
        'strand_2': '+'
    }}).T
    output_sv = normalize_sv_table(input_sv)
    assert (output_sv != expected_sv).sum().sum() == 0, output_sv

def test_pull_sv_supporting_reads_from_bundle(bundle):
    sv = {'chromosome_1':'chr1', 'position_1':1000, 'strand_1':'+',
          'chromosome_2':'chr2', 'position_2':2000, 'strand_2':'-',}
    sv_reads = pull_sv_supporting_reads_from_bundle(sv, bundle)
    expected_output = [bundle[0]]
    assert sv_reads == expected_output, (expected_output, sv_reads)

def test_find_presence_of_matching_sv():
    ix_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2']
    sv1 = pd.DataFrame([
        ['chr1', 1000, '+', 'chr2', 2000, '+'], 
        ['chr3', 3000, '+', 'chr4', 4000, '+'], 
    ], columns=ix_cols)
    sv2 = pd.DataFrame([
        ['chr1', 1049, '+', 'chr2', 2049, '+'], 
        ['chr3', 3000, '+', 'chr4', 4000, '-'], 
    ], columns=ix_cols)
    expected = pd.DataFrame([
        ['chr1', 1000, '+', 'chr2', 2000, '+', True], 
        ['chr3', 3000, '+', 'chr4', 4000, '+', False], 
    ], columns=ix_cols+['match'])
    sv1['match'] = find_presence_of_matching_sv(sv1, sv2, margin=50)
    assert (sv1 != expected).sum().sum() == 0, sv1

def test_extract_read_data_1():
    bam_path = 'tests/data/test.bam'
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1, end=2)
    assert df.shape[0] == 0, df

def test_extract_read_data_2():
    bam_path = 'tests/data/test.bam'
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1470, end=1477)
    assert df.shape[0] == 0, df

def test_extract_read_data_3():
    expected = np.array([['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'chr10', 51339923,
        51340887, '+', 35, 960, 3763, 35],
       ['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'PBEF1NeoTransposon',
        1478, 4996, '-', 288, 3480, 990, 990],
       ['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'chr10', 51340883,
        51341166, '+', 4466, 281, 11, 4466]])
    bam_path = 'tests/data/test.bam'
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1477, end=1478)
    assert np.all(df.to_numpy().astype(str) == expected.astype(str))
