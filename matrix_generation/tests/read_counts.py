import os 
import glob
import argparse
import sys
import gzip
import yaml

parent_path = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)  # gets the parent directory of the file that this script is in
sys.path.append(parent_path)  # adds the parent directory to the path

def read_counts_function_custom_folders(path_downsampled, path_full, output_directory):
    read_counts_full = {}
    file_list_full = glob.glob(f"{path_full}/*")
    baseline_file = f"{output_directory}/baseline.txt"
    
    print("starting baseline")
    
    if not os.path.exists(baseline_file):
        for file in file_list_full:
            base_name = os.path.basename(file)
            count = 0
            with gzip.open(file, "rt") as f:
                for _ in f:
                    count += 1
            read_counts_full[base_name] = int(count/4)

        with open(baseline_file, "w") as f:
            for k, v in read_counts_full.items():
                f.write(f"{k}: {v}\n")
    
    else:
        with open(f"{baseline_file}", "r") as f:
            for line in f:
                # Remove any leading and trailing whitespaces (like newline characters)
                line = line.strip()

                # Split the line by ': ' into key and value parts
                key, value = line.split(': ')

                # Add the key-value pair to the dictionary
                read_counts_full[key] = int(value)  # Convert the value to integer before storing

    print(read_counts_full)  # Print the dictionary to verify it's as expected
    
    print("done with baseline")

    read_counts_downsampled = {}
    
    file_list_downsampled = glob.glob(f"{path_downsampled}/*")
    for file in file_list_downsampled:
        base_name = os.path.basename(file)
        count = 0
        with gzip.open(file, "rt") as f:
            for _ in f:
                count += 1
        read_counts_downsampled[base_name] = int(count/4)    
    
    with open(f"{output_directory}/{os.path.basename(path_downsampled)}.txt", "w") as f:
        for k, v in read_counts_full.items():
            ratio_statement = f"{os.path.basename(path_downsampled)} {base_name} RATIO: {read_counts_downsampled[base_name] / read_counts_full[base_name]}"
            f.write(f"{ratio_statement}\n")
            print(f"{ratio_statement}")

if __name__ == "__main__":
    from fastq_processor import FastqProcessor
    
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    path_downsampled = fastq_processor.downsampled_fastq_directory
    path_full = fastq_processor.original_fastq_directory
    output_directory = f"{os.path.dirname(os.path.realpath(__file__))}/read_count_files_NEW"  # change ending for unique output

    
    if path_downsampled is not None or path_full is not None:
        with os.scandir(path_downsampled) as entries:
            for entry in entries:
                if entry.is_dir():
                    path_downsampled_specific = f"{path_downsampled}/{entry.name}"
                    read_counts_function_custom_folders(path_downsampled_specific, path_full, output_directory)
    print("done")