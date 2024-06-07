import numpy as np
import pandas as pd

from ontmont.process import merge_cn_segments, filter_sv_by_clone

def test_get_blocks_from_cn_1():
    cols = ['clone_id', 'chr', 'start', 'end', 'state']
    data = [
        ['A', '1', 10000, 20000, 2],
        ['A', '1', 30000, 40000, 2],
        ['A', '1', 40000, 50000, 1],
        ['A', '1', 50000, 60000, 2],
        ['A', '1', 260000, 300000, 2],
    ]
    cn = pd.DataFrame(columns=cols, data=data)
    merged = merge_cn_segments(cn, merge_gap=1000)
    expected = cn.copy()
    assert np.all(expected == merged)

def test_get_blocks_from_cn_2():
    cols = ['clone_id', 'chr', 'start', 'end', 'state']
    data = [
        ['A', '1', 10000, 20000, 2],
        ['A', '1', 30000, 40000, 2],
        ['A', '1', 40000, 50000, 1],
        ['A', '1', 50000, 60000, 2],
        ['A', '1', 260000, 300000, 2],
    ]
    cn = pd.DataFrame(columns=cols, data=data)
    merged = merge_cn_segments(cn, merge_gap=100000)
    expected = pd.DataFrame([
        ['A', '1', 10000, 40000, 2],
        ['A', '1', 40000, 50000, 1],
        ['A', '1', 50000, 60000, 2],
        ['A', '1', 260000, 300000, 2],
    ], columns=cols)
    assert np.all(expected == merged)

def test_get_blocks_from_cn_3():
    cols = ['clone_id', 'chr', 'start', 'end', 'state']
    data = [
        ['A', '1', 10000, 20000, 2],
        ['A', '2', 30000, 40000, 2],
        ['A', '2', 40000, 50000, 1],
        ['A', '2', 50000, 60000, 2],
        ['A', '2', 260000, 300000, 2],
    ]
    cn = pd.DataFrame(columns=cols, data=data)
    merged = merge_cn_segments(cn, merge_gap=500000)
    expected = pd.DataFrame([
        ['A', '1', 10000, 20000, 2],
        ['A', '2', 30000, 40000, 2],
        ['A', '2', 40000, 50000, 1],
        ['A', '2', 50000, 300000, 2],
    ], columns=cols)
    assert np.all(expected == merged)

def test_filter_sv_by_clone():
    cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 
            'vaf_A', 'vaf_B', 'vaf_C']
    clone_ids = ['A', 'B', 'C']
    data = [
        ['1', 1000, '+', '2', 2000, '-', 0.2, 0.0, 0.0],
        ['2', 1000, '+', '3', 2000, '-', 0.0, 0.3, 0.0],
        ['3', 1000, '+', '3', 2000, '-', 0.0, 0.0, 0.4],
        ['4', 1000, '+', '5', 4000, '-', 0.1, 0.0, 0.4],
    ]
    sv = pd.DataFrame(data=data, columns=cols)
    clonal_sv = filter_sv_by_clone(sv)
    expected = {
        'A': pd.DataFrame([['1', 1000, '+', '2', 2000, '-', 0.2, 0.0, 0.0]], columns=cols),
        'B': pd.DataFrame([['2', 1000, '+', '3', 2000, '-', 0.0, 0.3, 0.0]], columns=cols),
        'C': pd.DataFrame([['3', 1000, '+', '3', 2000, '-', 0.0, 0.0, 0.4]], columns=cols),
    }
    for clone_id in clone_ids:
        assert np.all(expected[clone_id] == clonal_sv[clone_id])