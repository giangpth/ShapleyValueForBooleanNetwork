import os
import re
import glob
import pandas as pd
import numpy as np
from sklearn.metrics import ndcg_score

def parse_log_file(filepath):
    """
    Parses a log file to extract predicted and ground-truth KO and KI values.
    Returns two dictionaries: pred_data and true_data.
    """
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Regex pattern to match lines like: NodeName       KO:     0.0 | KI:      -2.7
    pattern = re.compile(r'^([A-Za-z0-9_]+)\s+KO:\s+([-\d.]+)\s+\|\s+KI:\s+([-\d.]+)$', re.MULTILINE)
    matches = pattern.findall(content)

    numinputpattern = re.compile(r'NUMBER OF INPUT:\s+(\d+)')
    numinputmatch = numinputpattern.search(content)
    # print("numinputmatch: {}".format(numinputmatch))
    numinput = int(numinputmatch.group(1)) if numinputmatch else 0
    
    if len(matches) < 2:  # Expecting 2 blocks of 16 nodes
        return None, None

    nodes_per_block = len(matches) // 2
    
    pred_raw = matches[:nodes_per_block]
    true_raw = matches[nodes_per_block:2*nodes_per_block]

    pred_data = {node: {'KO': float(ko), 'KI': float(ki)} for node, ko, ki in pred_raw}
    true_data = {node: {'KO': float(ko), 'KI': float(ki)} for node, ko, ki in true_raw}

    return pred_data, true_data, numinput

def calculate_ndcg(true_dict, pred_dict, metric_key, numinput):
    """
    Filters nodes with non-zero values, handles negative values via shifting,
    handles single-item edge cases, and computes NDCG.
    """
    nodes = list(true_dict.keys())
    
    true_vals = np.array([true_dict[n][metric_key] for n in nodes])
    pred_vals = np.array([pred_dict[n][metric_key] for n in nodes])

    if len(set(true_vals)) == 1:
        # All true values are the same (e.g., all zeros), NDCG is undefined
        return 0.0, 0.0, -1.0
    
    # Filter nodes where at least one has a non-zero value
    non_zero_mask = (true_vals != 0.0) | (pred_vals != 0.0)

    nonzero_proportion = float(np.sum(non_zero_mask) / len(nodes))
    
    if not np.any(non_zero_mask):
        return 0.0, 0.0, -1.0
        
    y_true = (true_vals[non_zero_mask])
    y_pred = (pred_vals[non_zero_mask])

    rRMSE = float(np.sqrt(np.mean((y_true - y_pred) ** 2))/(2*(numinput+1)))
    # rRMSE = rRMSE/len(y_true)
    
    # If there is only 1 node remaining, NDCG cannot compute rankings (perfect match by default)
    if len(y_true) <= 1:
        if y_true[0] != 0 and y_pred[0] != 0:
            return 1.0, nonzero_proportion, rRMSE  # Perfect match, rRMSE is 0
        else:
            return 0.00000001, nonzero_proportion, rRMSE  # No match, rRMSE is -1

    # Shift values to be non-negative if there are negative numbers (scikit-learn requirement)
    min_val = min(y_true.min(), y_pred.min())
    if min_val < 0:
        shift = abs(min_val)
        y_true = y_true + shift
        y_pred = y_pred + shift
        
    try:
        score = ndcg_score([y_true], [y_pred])
        if score == 0.0:
            score = 0.00000001  # Avoid exact zero for downstream calculations
    except ValueError:
        print(f"Error computing NDCG for true: {y_true}, pred: {y_pred}")
        print(f"Error computing rRMSE for true: {y_true}, pred: {y_pred}")
        if len(y_true) == 1 and len(y_pred) == 1:
            score = 1.0 if y_true[0] != 0 and y_pred[0] != 0 else 0.00000001
        else:
            print(f"Error computing NDCG for true: {y_true}, pred: {y_pred}")
            score = np.nan

        
    return score, nonzero_proportion, rRMSE

def compute_ndcg_statistics_across_runs(folder_paths):
    """
    Iterates over 5 run folders, calculates NDCG for each file per run,
    then computes the mean and variance of the NDCG scores across the 5 runs for each log file.
    """
    # Dictionary to store NDCG scores per file across runs:
    # ndcg_records[filename] = {'KO': [run1_score, run2_score, run3_score, run4_score, run5_score], 'KI': [...]}
    ndcg_records = {}

    for folder in folder_paths:
        if not os.path.exists(folder):
            print(f"Warning: Folder not found: {folder}")
            continue
        else:
            print("Processing folder {}".format(folder))
            
        log_files = glob.glob(os.path.join(folder, "*.log"))
        
        for filepath in log_files:
            filename = os.path.basename(filepath)
            print("Processing file {}".format(filename))
            pred_data, true_data, numinput = parse_log_file(filepath)
            
            if pred_data is None or true_data is None:
                continue

            print(calculate_ndcg(true_data, pred_data, 'KO', numinput))
            print(calculate_ndcg(true_data, pred_data, 'KI', numinput))
            ndcg_ko, nonzero_ko, rRMSE_ko = calculate_ndcg(true_data, pred_data, 'KO', numinput)
            ndcg_ki, nonzero_ki, rRMSE_ki = calculate_ndcg(true_data, pred_data, 'KI', numinput)

            if filename not in ndcg_records:
                ndcg_records[filename] = {'KO': [], 'KI': []}
                ndcg_records[filename]['KO_nonzero'] = []
                ndcg_records[filename]['KI_nonzero'] = []
                ndcg_records[filename]['KO_rRMSE'] = []
                ndcg_records[filename]['KI_rRMSE'] = []
                
            ndcg_records[filename]['KO'].append(ndcg_ko)
            ndcg_records[filename]['KI'].append(ndcg_ki)
            ndcg_records[filename]['KO_nonzero'].append(nonzero_ko)
            ndcg_records[filename]['KI_nonzero'].append(nonzero_ki)
            ndcg_records[filename]['KO_rRMSE'].append(rRMSE_ko)
            ndcg_records[filename]['KI_rRMSE'].append(rRMSE_ki)

    # Compute mean and variance across the runs for each file
    summary_list = []

    for filename, scores in ndcg_records.items():
        ko_scores = scores['KO']
        ki_scores = scores['KI']
        rMSE_ko_scores = scores['KO_rRMSE']
        print(rMSE_ko_scores)
        rMSE_ki_scores = scores['KI_rRMSE']
        ko_nonzero = scores['KO_nonzero']
        ki_nonzero = scores['KI_nonzero']
        num_runs = len(ko_scores)
        
        summary_list.append({
            'Log_File': filename,
            'Num_Runs': num_runs,
            
            # KO NDCG Stats
            'NDCG_KO_Mean': np.mean(ko_scores),
            'NDCG_KO_Var': np.var(ko_scores, ddof=1) if num_runs > 1 else 0.0,
            
            # KI NDCG Stats
            'NDCG_KI_Mean': np.mean(ki_scores),
            'NDCG_KI_Var': np.var(ki_scores, ddof=1) if num_runs > 1 else 0.0,

            'Nonzero_KO_mean': np.mean(ko_nonzero),
            'Nonzero_KO_Var': np.var(ko_nonzero, ddof=1) if num_runs > 1 else 0.0,

            'Nonzero_KI_mean': np.mean(ki_nonzero),
            'Nonzero_KI_Var': np.var(ki_nonzero, ddof=1) if num_runs > 1 else 0.0,


            'rRMSE_KO': np.mean(rMSE_ko_scores),
            'rRMSE_KO_Var': np.var(rMSE_ko_scores, ddof=1) if num_runs > 1 else 0.0,

            'rRMSE_KI': np.mean(rMSE_ki_scores),
            'rRMSE_KI_Var': np.var(rMSE_ki_scores, ddof=1) if num_runs > 1 else 0.0
        })

    return pd.DataFrame(summary_list)

if __name__ == "__main__":
    # Specify paths to your 3 run folders here
    folders = ['logs/realbenchmark_propagation_1strun', 'logs/realbenchmark_propagation_2ndRun', 'logs/realbenchmark_propagation_3rdRun', 'logs/realbenchmark_propagation_4thRun', 'logs/realbenchmark_propagation_5thRun']
    
    print("Calculating NDCG mean and variance across 5 runs...")
    df_ndcg_stats = compute_ndcg_statistics_across_runs(folders)
    
    if not df_ndcg_stats.empty:
        print("\n--- NDCG Statistics Across 5 Runs ---")
        print(df_ndcg_stats.to_string(index=False))
        
        # Save summary to CSV
        output_csv = 'ndcg_stats_across_5_runs.csv'
        df_ndcg_stats.to_csv(output_csv, index=False)
        print(f"\nResults successfully saved to '{output_csv}'")
    else:
        print("No valid log files processed.")