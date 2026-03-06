import sys
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from patsy import dmatrix
from scipy.interpolate import UnivariateSpline

data_pol2 = pd.read_csv('./predicted_rN_pol2.tsv', header = None, sep="\t")
data_pol3 = pd.read_csv('./predicted_rN_pol3.tsv', header = None, sep="\t")
rep_ori = pd.read_csv('/rep_ori.tsv', header = None, sep="\t")
rep_ori.columns = ['status', 'chr', 'start', 'end', 'contig']

data_pol2.columns = [
    'read_name',
    'contig',
    'position',
    'strand',
    'ampl_mean_2nd',
    'ampl_mean_3rd',
    'ampl_mean_4th',
    'ampl_mean_5th',
    'ampl_mean_6th',
    'ampl_mean_7th',
    'model_deviation_2nd',
    'model_deviation_3rd',
    'model_deviation_4th',
    'model_deviation_5th',
    'model_deviation_6th',
    'model_deviation_7th',
    'rType_pred',
    'rType_prob',
    'base_consensus'
]

data_pol3.columns = [
    'read_name',
    'contig',
    'position',
    'strand',
    'ampl_mean_2nd',
    'ampl_mean_3rd',
    'ampl_mean_4th',
    'ampl_mean_5th',
    'ampl_mean_6th',
    'ampl_mean_7th',
    'model_deviation_2nd',
    'model_deviation_3rd',
    'model_deviation_4th',
    'model_deviation_5th',
    'model_deviation_6th',
    'model_deviation_7th',
    'rType_pred',
    'rType_prob',
    'base_consensus'
]

data_pol2['rType_prob'] = pd.to_numeric(data_pol2['rType_prob'], errors='coerce')
data_pol3['rType_prob'] = pd.to_numeric(data_pol3['rType_prob'], errors='coerce')

data_pol2_filtered = data_pol2[(data_pol2['rType_prob']>=0.99)]
data_pol3_filtered = data_pol3[(data_pol3['rType_prob']>=0.99)]

rType_mapping = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
rType_mapping_rev = {1: 'T', 2: 'G', 3: 'C', 4: 'A'}
data_pol2_filtered.loc[data_pol2_filtered['strand']=="pos", 'rType_pred_base'] = data_pol2_filtered.loc[data_pol2_filtered['strand']=="pos", 'rType_pred'].map(rType_mapping)
data_pol2_filtered.loc[data_pol2_filtered['strand']=="rev", 'rType_pred_base'] = data_pol2_filtered.loc[data_pol2_filtered['strand']=="rev", 'rType_pred'].map(rType_mapping_rev)
data_pol3_filtered.loc[data_pol3_filtered['strand']=="pos", 'rType_pred_base'] = data_pol3_filtered.loc[data_pol3_filtered['strand']=="pos", 'rType_pred'].map(rType_mapping)
data_pol3_filtered.loc[data_pol3_filtered['strand']=="rev", 'rType_pred_base'] = data_pol3_filtered.loc[data_pol3_filtered['strand']=="rev", 'rType_pred'].map(rType_mapping_rev)
data_pol2_filtered = data_pol2_filtered[(data_pol2_filtered['rType_pred_base']==data_pol2_filtered['base_consensus'])]
data_pol3_filtered = data_pol3_filtered[(data_pol3_filtered['rType_pred_base']==data_pol3_filtered['base_consensus'])]


grouped_counts_pol2 = data_pol2_filtered.groupby(['contig', 'position', 'strand']).size().reset_index(name='count')
grouped_counts_pol3 = data_pol3_filtered.groupby(['contig', 'position', 'strand']).size().reset_index(name='count')

grouped_counts_pol2['position'] = pd.to_numeric(grouped_counts_pol2['position'], errors='raise')
grouped_counts_pol3['position'] = pd.to_numeric(grouped_counts_pol3['position'], errors='raise')

pivoted_pol2 = grouped_counts_pol2.pivot(index=['contig', 'position'],columns = 'strand', values='count')
pivoted_pol2 = pivoted_pol2.rename(columns={'+': 'count_pos', '-': 'count_rev'}).reset_index()
pivoted_pol2 = pivoted_pol2.fillna(0)

pivoted_pol3 = grouped_counts_pol2.pivot(index=['contig', 'position'],columns = 'strand', values='count')
pivoted_pol3 = pivoted_pol3.rename(columns={'+': 'count_pos', '-': 'count_rev'}).reset_index()
pivoted_pol3 = pivoted_pol3.fillna(0)

def add_rolling_sums_per_contig(df, half_window=1000):
    # df is one contig subset
    df = df.sort_values('position').reset_index(drop=True)

    positions = df['position'].to_numpy()
    n = len(df)

    rolling_sum_pos = np.empty(n, dtype=float)
    rolling_sum_rev = np.empty(n, dtype=float)

    start_idx = 0
    end_idx = 0

    pos_arr = df['pos'].to_numpy()
    rev_arr = df['rev'].to_numpy()

    # (optional) prefix sums make each window sum O(1), overall O(n)
    pos_cum = np.concatenate(([0.0], np.cumsum(pos_arr)))
    rev_cum = np.concatenate(([0.0], np.cumsum(rev_arr)))

    for i in range(n):
        center = positions[i]
        window_start = center - half_window
        window_end = center + half_window

        while start_idx < n and positions[start_idx] < window_start:
            start_idx += 1
        while end_idx < n and positions[end_idx] <= window_end:
            end_idx += 1

        rolling_sum_pos[i] = pos_cum[end_idx] - pos_cum[start_idx]
        rolling_sum_rev[i] = rev_cum[end_idx] - rev_cum[start_idx]

    df['rolling_sum_pos'] = rolling_sum_pos
    df['rolling_sum_rev'] = rolling_sum_rev

    denom = df['rolling_sum_pos'] + df['rolling_sum_rev']
    df['rolling_median_enrichment'] = np.where(
        denom > 0,
        (df['rolling_sum_pos'] - df['rolling_sum_rev']) / denom,
        np.nan
    )

    return df


pivoted_pol2 = (
    pivoted_pol2
    .groupby('contig', group_keys=False)
    .apply(add_rolling_sums_per_contig, half_window=1000)
)

pivoted_pol3 = (
    pivoted_pol3
    .groupby('contig', group_keys=False)
    .apply(add_rolling_sums_per_contig, half_window=1000)
)

min_gap = 20000

confirmed = (
    rep_ori.loc[rep_ori['status'].isin(["Confirmed","Likely"]), ['contig', 'start', 'end']]
    .dropna(subset=['contig', 'start', 'end'])
    .assign(start=lambda d: d['start'].astype(int),
            end=lambda d: d['end'].astype(int))
    .sort_values(['contig', 'start'])
)

def merge_close_intervals(g, min_gap=5000):
    starts = g['start'].to_numpy()
    ends   = g['end'].to_numpy()

    merged_starts = []
    merged_ends = []

    cur_start = starts[0]
    cur_end = ends[0]

    for s, e in zip(starts[1:], ends[1:]):
        if s <= cur_end + min_gap:
            cur_end = max(cur_end, e)
        else:
            merged_starts.append(cur_start)
            merged_ends.append(cur_end)
            cur_start, cur_end = s, e

    merged_starts.append(cur_start)
    merged_ends.append(cur_end)

    return pd.DataFrame({
        'contig': g['contig'].iloc[0],
        'start': merged_starts,
        'end': merged_ends
    })


def drop_close_pairs(g):
    next_start = g['start'].shift(-1)
    gap_to_next = next_start - g['end']  

    bad_pair = gap_to_next.le(min_gap).fillna(False)

    drop = bad_pair | bad_pair.shift(1, fill_value=False)

    return g.loc[~drop]
    
def drop_edges(g, n=1):
    if len(g) <= 2 * n:
        return g.iloc[0:0]
    return g.iloc[n:-n]

confirmed_filtered = (
    confirmed
    .groupby('contig', group_keys=False)
    .apply(merge_close_intervals, min_gap=min_gap) 
    .reset_index(drop=True) 
    .sort_values(['contig','start']) 
    .groupby('contig', group_keys=False)
    .apply(drop_edges, n=2) 
)

pivoted_pol2 = pivoted_pol2.sort_values(['contig', 'position']).reset_index(drop=True)
pivoted_pol3 = pivoted_pol3.sort_values(['contig', 'position']).reset_index(drop=True)

contig_arrays = {}
for contig, g in pivoted_pol2.groupby('contig', sort=False):
    pos = g['position'].to_numpy()
    enr = g['rolling_median_enrichment'].to_numpy()
    contig_arrays[contig] = (pos, enr)

all_chunks = []

W = 20000

for contig, ev in confirmed_filtered.groupby('contig', sort=False):
    if contig not in contig_arrays:
        continue

    pos, enr = contig_arrays[contig]

    for start, end in ev[['start', 'end']].itertuples(index=False, name=None):

        l = np.searchsorted(pos, start - W, side='left')
        r = np.searchsorted(pos, start,     side='right')
        if r > l:
            rp = pos[l:r] - start
            all_chunks.append(pd.DataFrame({
                'relative_pos': rp,
                'rolling_median_enrichment': enr[l:r]
            }))

        l = np.searchsorted(pos, end,     side='left')
        r = np.searchsorted(pos, end + W, side='right')
        if r > l:
            rp = pos[l:r] - end
            all_chunks.append(pd.DataFrame({
                'relative_pos': rp,
                'rolling_median_enrichment': enr[l:r]
            }))

stacked = pd.concat(all_chunks, ignore_index=True)

agg_pol2 = (
    stacked.groupby('relative_pos', sort=True)['rolling_median_enrichment']
    .agg(['mean', 'median'])
    .reset_index()
)

agg_filtered_pol2 = agg_pol2.query("-5000 <= relative_pos <= 5000")

contig_arrays = {}
for contig, g in pivoted_pol3.groupby('contig', sort=False):
    pos = g['position'].to_numpy()
    enr = g['rolling_median_enrichment'].to_numpy()
    contig_arrays[contig] = (pos, enr)

all_chunks = []

W = 20000

for contig, ev in confirmed_filtered.groupby('contig', sort=False):
    if contig not in contig_arrays:
        continue

    pos, enr = contig_arrays[contig]

    for start, end in ev[['start', 'end']].itertuples(index=False, name=None):

        l = np.searchsorted(pos, start - W, side='left')
        r = np.searchsorted(pos, start,     side='right')
        if r > l:
            rp = pos[l:r] - start
            all_chunks.append(pd.DataFrame({
                'relative_pos': rp,
                'rolling_median_enrichment': enr[l:r]
            }))

        l = np.searchsorted(pos, end,     side='left')
        r = np.searchsorted(pos, end + W, side='right')
        if r > l:
            rp = pos[l:r] - end
            all_chunks.append(pd.DataFrame({
                'relative_pos': rp,
                'rolling_median_enrichment': enr[l:r]
            }))

stacked = pd.concat(all_chunks, ignore_index=True)

agg_pol3 = (
    stacked.groupby('relative_pos', sort=True)['rolling_median_enrichment']
    .agg(['mean', 'median'])
    .reset_index()
)

agg_filtered_pol3 = agg_pol3.query("-5000 <= relative_pos <= 5000")

agg_filtered_pol2 = agg_filtered_pol2[
    agg_filtered_pol2['mean'].abs() <= 0.5 #Remove outliers
]
agg_filtered_pol3 = agg_filtered_pol3[
    agg_filtered_pol3['mean'].abs() <= 0.5 #Remove outliers
]
plt.figure(figsize=(5, 5))
x2 = agg_filtered_pol2['relative_pos'].values
y2 = agg_filtered_pol2['mean'].values
spline_pol2 = UnivariateSpline(x2, y2, s=15.8)
plt.plot(x2, spline_pol2(x2), color='#698a9c', linestyle='-', linewidth=2, label='pol2-M644G')

x3 = agg_filtered_pol3['relative_pos'].values
y3 = agg_filtered_pol3['mean'].values
spline_pol3 = UnivariateSpline(x3, y3, s=15.2)
plt.plot(x3, spline_pol3(x3), color='#d47277', linestyle='-', linewidth=2, label='pol3-L612M')

plt.axvline(0, color='red', linestyle='--', linewidth=1)
plt.xlabel("Relative Position (bp)")
plt.ylabel("rNMP Enrichment on Watson Strand")
plt.xticks([-5000, -3000, -1000, 0, 1000, 3000, 5000])
plt.yticks([-0.05, -0.04, -0.03, -0.02, -0.01, 
             0,
             0.01, 0.02, 0.03, 0.04, 0.05])
plt.legend()
plt.ylim(-0.1, 0.1)
plt.savefig("../model_training/data/ygDNA_ribo_enrichment_XGB.pdf")