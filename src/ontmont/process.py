import pandas as pd

def update_blocks_and_reset_prev(blocks, prev, row, 
        features=['clone_id', 'chr', 'start', 'end', 'state']):
    if prev['chr']: # not null
        blocks = pd.concat([blocks, prev.to_frame().T])
    prev = pd.Series({f: row[f] for f in features}, index=features)
    return blocks, prev
    
def merge_cn_segments(cn, merge_gap=100000):
    #in_pp = '/juno/work/shah/isabl_data_lake/analyses/37/74/23774/SPECTRUM-OV-081_S1_LEFT_OVARY_cn.csv'
    features = ['clone_id', 'chr', 'start', 'end', 'state']
    cn = cn[features].copy() # reduce features for the sake of memory

    # merge segments with same metric
    prev = pd.Series({f:None for f in features}, index=features) # init prev cn
    blocks = pd.DataFrame(columns=features) # init empty blocks

    for rix, row in cn.iterrows():
        if prev['chr'] != row['chr']:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if (prev['state'] != row['state']):
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        if (prev['state'] == row['state']) and abs(prev['end'] - row['start']) > merge_gap:
            blocks, prev = update_blocks_and_reset_prev(blocks, prev, row)
            continue
        # update prev block
        prev['end'] = row['end']

    blocks = pd.concat([blocks, prev.to_frame().T]) # add prev
    blocks.loc[:, ['start', 'end', 'state']] = blocks.loc[:, ['start', 'end', 'state']].astype(int)
    blocks['start'] = blocks['start'] # so segments don't overlap; w/o this !mixFlag[j] error occurs in MutationTimeR step
    return blocks.reset_index(drop=True)

def filter_sv_by_clone(sv, vaf_col_str='vaf_'):
    vaf_cols = [c for c in sv.columns if c.startswith(vaf_col_str)]
    clone_ids = [c.replace(vaf_col_str, '') for c in vaf_cols]
    clonal_sv = {}
    for clone_id in clone_ids:
        clonal = f'vaf_{clone_id}'
        others = [f'vaf_{c}' for c in clone_ids if c != clone_id]
        clonal_sv[clone_id] = sv[
            (sv[clonal] > 0) &
            (sv[others].sum(axis=1) == 0)
        ].reset_index(drop=True)
    return clonal_sv