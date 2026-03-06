import cudf
import joblib
import numpy as np
import pandas as pd
from xgboost import XGBClassifier
import sys
import os
import glob
from tqdm import tqdm

for d in ["numba_cache", "cupy_cache", "xdg_cache", "tmp"]:
    os.makedirs(d, exist_ok=True)
os.environ['NUMBA_CACHE_DIR'] = './numba_cache'
os.environ['CUPY_CACHE_DIR'] = './cupy_cache'
os.environ['XDG_CACHE_HOME'] = './xdg_cache'
os.environ['TMPDIR'] = './tmp'

chunk_path = sys.argv[1]  
rn_for_path = sys.argv[2]
rn_rev_path = sys.argv[3]
model_path = sys.argv[4]
coord_map_path = sys.argv[5]
output_path = sys.argv[6]

# Read in strand information for each read_name
rn_for = cudf.read_csv(rn_for_path, sep='\t', header = None)
rn_rev = cudf.read_csv(rn_rev_path, sep='\t', header = None)
rn_for.columns = ['read_name', 'strand']
rn_rev.columns = ['read_name', 'strand']
rn_strand = cudf.concat([rn_for, rn_rev], ignore_index=True)

# Read in the trained model
final_model = joblib.load(model_path)
features = ["ampl_mean_2nd","ampl_mean_3rd","ampl_mean_4th", "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th", "model_deviation_2nd","model_deviation_3rd","model_deviation_4th", "model_deviation_5th", "model_deviation_6th", "model_deviation_7th"]

# Read in the original coordinates of nucleotides on the official ref gn
origi_coord = cudf.read_csv(coord_map_path, sep = " ", header = None)
origi_coord.columns = ['contig', 'base', 'position']

# Read in data in chunks and start preparation and prediction
chunk_files = glob.glob(os.path.join('./', chunk_path, 'chunk_*.txt'))
for file in tqdm(chunk_files, desc="Processing files", unit="file"):
    #print(f"Reading from: {file}")


    # Read in eventalign data
    eventalign_data = cudf.read_csv(file, sep='\t', usecols=["contig", "position", "reference_kmer", "read_name", "event_level_mean", "event_length", "model_kmer", "model_mean"])

    # Filter out rows where model_kmer is "NNNNNNNNN"
    eventalign_mut = eventalign_data[eventalign_data["model_kmer"] != "NNNNNNNNN"].copy()

    # Addmodel_deviation = event_level_mean - model_mean
    eventalign_mut["model_deviation"] = (
        eventalign_mut["event_level_mean"] - eventalign_mut["model_mean"]
    )

    # Calculate min/max positions per read_name
    position_bounds = (
        eventalign_mut
        .groupby("read_name")
        .agg({"position": ["min", "max"]})
    )

    # Flatten multi-index columns
    position_bounds.columns = ["min_pos", "max_pos"]
    position_bounds = position_bounds.reset_index()

    # Join bounds back to main dataframe
    eventalign_with_bounds = eventalign_mut.merge(position_bounds, on="read_name")

    # Filter out positions within 20 bases of the ends
    eventalign_trimmed = eventalign_with_bounds[
        (eventalign_with_bounds["position"] >= eventalign_with_bounds["min_pos"] + 20) &
        (eventalign_with_bounds["position"] <= eventalign_with_bounds["max_pos"] - 20)
    ]
    # Summarize by (position, read_name)
    eventalign_summary = (
        eventalign_trimmed
        .groupby(["contig","position", "read_name"])
        .agg({
            "event_level_mean": "mean",
            "model_deviation": "mean",
            "event_length": "sum"
        })
        .rename(columns={
            "event_level_mean": "event_level_mean_mean",
            "model_deviation": "model_deviation_mean",
            "event_length": "dwell_time"
        })
        .reset_index()
    )

    eventalign_summary_withstrand = eventalign_summary.merge(rn_strand, on = "read_name", how = "left")

    # Define kmer-nucleotide position mapping for strand == pos
    pos_shifts = list(zip(range(7, 1, -1), ["2nd", "3rd", "4th", "5th", "6th", "7th"]))

    # Define kmer-nucleotide position mapping for strand == rev
    rev_shifts = list(zip(range(6, 0, -1), ["7th", "6th", "5th", "4th", "3rd", "2nd"]))

    # Get k-mer data for forward reads
    df_pos = eventalign_summary_withstrand[eventalign_summary_withstrand["strand"] == "pos"]
    df_pos_ori = df_pos.copy()

    ## Filter valid range (start+7 to end-6 for each read to get all 2nd-7th kmer of a position)
    pos_bounds = df_pos.groupby("read_name").agg({"position": ["min", "max"]}).reset_index()
    pos_bounds.columns = ["read_name", "min_pos", "max_pos"]
    df_pos = df_pos.merge(pos_bounds, on="read_name")
    df_pos = df_pos[(df_pos["position"] >= df_pos["min_pos"] + 7) & (df_pos["position"] <= df_pos["max_pos"] - 6)]
    df_pos = df_pos.drop(columns=["min_pos", "max_pos"])

    ## Obtain kmer data
    result_pos = df_pos[["read_name","contig", "position", "strand"]].copy()

    for shift, label in pos_shifts:
        lag_df = df_pos_ori[["read_name", "contig", "position", "event_level_mean_mean", "model_deviation_mean"]].copy()
        lag_df["position"] = lag_df["position"] + shift  # shift backward to p

        lag_df = lag_df.rename(columns={
            "event_level_mean_mean": f"ampl_mean_{label}",
            "model_deviation_mean": f"model_deviation_{label}"
        })

        result_pos = result_pos.merge(lag_df, on=["read_name","contig", "position"], how="left")

    # Get k-mer data for reverse reads
    df_rev = eventalign_summary_withstrand[eventalign_summary_withstrand["strand"] == "rev"]
    df_rev_ori = df_rev.copy()

    ## Filter valid range (start+6 to end-7)
    rev_bounds = df_rev.groupby("read_name").agg({"position": ["min", "max"]}).reset_index()
    rev_bounds.columns = ["read_name", "min_pos", "max_pos"]
    df_rev = df_rev.merge(rev_bounds, on="read_name")
    df_rev = df_rev[
        (df_rev["position"] >= df_rev["min_pos"] + 6) &
        (df_rev["position"] <= df_rev["max_pos"] - 7)
    ].drop(columns=["min_pos", "max_pos"])

    ## Obtain kmer data
    result_rev = df_rev[["read_name","contig", "position", "strand"]].copy()

    for shift, label in rev_shifts:
        lag_df = df_rev_ori[["read_name","contig", "position", "event_level_mean_mean", "model_deviation_mean"]].copy()
        lag_df["position"] = lag_df["position"] + shift  # shift data from future (p + shift) to current p

        # Rename to match biological ordering
        lag_df = lag_df.rename(columns={
            "event_level_mean_mean": f"ampl_mean_{label}",
            "model_deviation_mean": f"model_deviation_{label}"
        })

        result_rev = result_rev.merge(lag_df, on=["read_name","contig", "position"], how="left")

    # Merge forward and reverse
    final_result = cudf.concat([result_pos, result_rev], ignore_index=True)

    ordered_columns = [
        "read_name", "contig", "position", "strand",
        "ampl_mean_2nd", "ampl_mean_3rd", "ampl_mean_4th",
        "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th",
        "model_deviation_2nd", "model_deviation_3rd", "model_deviation_4th",
        "model_deviation_5th", "model_deviation_6th", "model_deviation_7th"
    ]

    final_result = final_result[ordered_columns]
    final_result = final_result.sort_values(by=["read_name","contig", "position"]).reset_index(drop=True)

    # Do XGB prediction
    ## Prepare input
    final_model.set_params(device="cuda")
    final_result_cp = final_result.copy()
    final_result_cp[features] = final_result_cp[features].fillna(np.nan)
    final_result_cp[features] = final_result_cp[features].fillna(np.nan)
    X = final_result_cp[features].values.astype(np.float32)
    y_pred = final_model.predict(X)
    y_prob = final_model.predict_proba(X)
    y_pred_prob = y_prob[np.arange(len(y_pred)), y_pred]
    final_result['rType_pred'] = y_pred
    final_result['rType_prob'] = y_pred_prob
    final_result[features] = final_result[features].fillna(np.nan)

    # Get the coordinate on the official ref genome
    final_result_basein = final_result.merge(origi_coord, on = ["contig", "position"], how = 'left')
    final_result_basein.columns = [
        "read_name","contig", "position", "strand",
        "ampl_mean_2nd", "ampl_mean_3rd", "ampl_mean_4th",
        "ampl_mean_5th", "ampl_mean_6th", "ampl_mean_7th",
        "model_deviation_2nd", "model_deviation_3rd", "model_deviation_4th",
        "model_deviation_5th", "model_deviation_6th", "model_deviation_7th",
        "rType_pred", "rType_prob","base_consensus"
    ]


    # Save the result
    filename = os.path.basename(file)
    output = os.path.join('./', output_path, filename)
    final_result_basein.to_csv(output, sep='\t', index=False)



