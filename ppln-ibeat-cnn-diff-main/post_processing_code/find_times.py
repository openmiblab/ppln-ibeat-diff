import re
import pandas as pd
from datetime import datetime
import os

def analyze_test_dataset_performance(log_path):
    if not os.path.exists(log_path):
        print(f"Error: Could not find '{log_path}'")
        return

    data = []
    last_high_res_ts = None
    is_testing_phase = False

    with open(log_path, 'r') as f:
        for line in f:
            # Only start recording once we hit the test section of the log
            if "Testing Results" in line or "inference" in line.lower():
                is_testing_phase = True
            
            # Capture high-res timestamp
            ts_match = re.search(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})', line)
            if ts_match:
                last_high_res_ts = ts_match.group(1)
            
            # Pair timestamp with the specific saved slice
            if "Saved results for test case" in line:
                case_match = re.search(r'test case ([\w-]+)', line)
                if case_match and last_high_res_ts:
                    case_name = case_match.group(1)
                    data.append({'timestamp': last_high_res_ts, 'case': case_name})

    if not data:
        print("No test data found. Ensure the log includes the final inference/testing phase.")
        return

    df = pd.DataFrame(data)
    df['dt'] = pd.to_datetime(df['timestamp'], format='%Y-%m-%d %H:%M:%S,%f')
    df = df.sort_values('dt')
    df['time_diff'] = df['dt'].diff().dt.total_seconds()

    # Filter out gaps > 60s (breaks between datasets)
    inference_df = df[df['time_diff'] < 60].copy()

    # IMPROVED GROUPING: This looks for the unique Scan ID (e.g., iBE-2128-baseline)
    # to make sure we get 9 unique datasets even if the Patient ID is the same.
    def get_dataset_id(case):
        parts = case.split('_')
        # Join Patient ID and Scan Type (e.g., iBE-2128 + baseline)
        return "_".join(parts[:2])

    inference_df['dataset_id'] = inference_df['case'].apply(get_dataset_id)

    # Calculate Total Time per Full Dataset
    dataset_totals = inference_df.groupby('dataset_id')['time_diff'].sum().reset_index()
    dataset_totals.columns = ['Dataset ID', 'Total Time (s)']
    
    # Calculate Grand Mean
    mean_dataset_time = dataset_totals['Total Time (s)'].mean()

    print("\n" + "="*60)
    print(f"   PERFORMANCE FOR {len(dataset_totals)} WARPED TEST DATASETS")
    print("="*60)
    print(f"{'Dataset ID':<35} | {'Total Time (s)':<15}")
    print("-"*60)
    for _, row in dataset_totals.iterrows():
        print(f"{row['Dataset ID']:<35} | {row['Total Time (s)']:<15.2f}")
    
    print("-" * 60)
    print(f"MEAN TIME PER FULL DATASET: {mean_dataset_time:.2f} seconds")
    print(f"TOTAL CLINICAL THROUGHPUT:  {len(dataset_totals) / (dataset_totals['Total Time (s)'].sum() / 60):.2f} datasets/min")
    print("="*60)

if __name__ == "__main__":
    downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
    log_file = os.path.join(downloads_path, "train_with_times.log")
    analyze_test_dataset_performance(log_file)