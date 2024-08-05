import os
import pandas as pd
import argparse

def read_benchmark_file(filepath):
    """
    Reads a benchmark file and returns its content as a DataFrame.

    Args:
        filepath (str): Path to the benchmark file.

    Returns:
        pd.DataFrame: DataFrame with benchmark data.
    """
    return pd.read_csv(filepath, sep='\t')

def convert_time_to_hms(seconds):
    """
    Converts time from seconds to hh:mm:ss format.

    Args:
        seconds (float): Time in seconds.

    Returns:
        str: Time in hh:mm:ss format.
    """
    seconds = float(seconds)
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = int(seconds % 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"

def convert_memory_to_gb(mb):
    """
    Converts memory from MB to GB, rounds to 3 decimal places.
    If the value is below 1 GB, it will still be correctly formatted.

    Args:
        mb (float): Memory in MB.

    Returns:
        float: Memory in GB.
    """
    gb = float(mb) / 1024
    return round(gb, 3) if gb >= 1 else round(float(mb) / 1024, 3)

def collect_benchmarks(root_dir):
    """
    Collects benchmark data from the directory tree and returns a dictionary of DataFrames.

    Args:
        root_dir (str): Root directory containing benchmark files.

    Returns:
        dict: Dictionary with rule names as keys and DataFrames as values.
    """
    data = {}

    for subject_dir in os.listdir(root_dir):
        subject_path = os.path.join(root_dir, subject_dir)
        if os.path.isdir(subject_path):
            for benchmark_file in os.listdir(subject_path):
                if benchmark_file.endswith('.benchmark.txt'):
                    rule = benchmark_file.replace('.benchmark.txt', '')
                    filepath = os.path.join(subject_path, benchmark_file)
                    benchmark_data = read_benchmark_file(filepath)
                    benchmark_data['subject'] = subject_dir
                    if rule not in data:
                        data[rule] = []
                    data[rule].append(benchmark_data)

    # Combine lists of DataFrames into single DataFrames per rule
    rule_dfs = {rule: pd.concat(df_list, ignore_index=True) for rule, df_list in data.items()}
    return rule_dfs

def create_summary_table(rule_dfs):
    """
    Creates a summary table with rules as rows and benchmark modalities as columns.
    Converts time to hh:mm:ss and memory to GB.

    Args:
        rule_dfs (dict): Dictionary of DataFrames with rule names as keys.

    Returns:
        pd.DataFrame: Summary table.
    """
    summary_data = []

    for rule, df in rule_dfs.items():
        rule_summary = {
            'rule': rule,
            'max_rss (GB)': convert_memory_to_gb(df['max_rss'].mean()),
            'max_vms (GB)': convert_memory_to_gb(df['max_vms'].mean()),
            'max_uss (GB)': convert_memory_to_gb(df['max_uss'].mean()),
            'max_pss (GB)': convert_memory_to_gb(df['max_pss'].mean()),
            'io_in': round(df['io_in'].mean(), 3),
            'io_out': round(df['io_out'].mean(), 3),
            'mean_load': round(df['mean_load'].mean(), 3),
            'cpu_time (hh:mm:ss)': convert_time_to_hms(df['cpu_time'].mean())
        }
        summary_data.append(rule_summary)

    summary_df = pd.DataFrame(summary_data)
    summary_df.set_index('rule', inplace=True)
    return summary_df

def main():
    """
    Main function to execute the script.
    """
    parser = argparse.ArgumentParser(description='Process benchmark files.')
    parser.add_argument('benchmark_dir', type=str, help='Root directory containing benchmark files.')
    args = parser.parse_args()
    
    root_dir = args.benchmark_dir
    rule_dfs = collect_benchmarks(root_dir)

    for rule, df in rule_dfs.items():
        # Ensure subject is the first column
        cols = ['subject'] + [col for col in df.columns if col != 'subject']
        df = df[cols]
        output_file = os.path.join(root_dir, f'{rule}_benchmark_results.csv')
        df.to_csv(output_file, index=False)
        print(f"Benchmark results for {rule} have been saved to {output_file}")

    summary_df = create_summary_table(rule_dfs)
    summary_output_file = os.path.join(root_dir, 'summary_benchmark_results.csv')
    summary_df.to_csv(summary_output_file)
    print(f"Summary benchmark results have been saved to {summary_output_file}")

if __name__ == "__main__":
    main()