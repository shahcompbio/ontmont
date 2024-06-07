import pytest

import yaml
import pyfaidx
import pandas as pd
import numpy as np

from ontmont.utils import remove_duplicates_from_tra_table, filter_sv_with_breakpoint_at_contig_ends


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