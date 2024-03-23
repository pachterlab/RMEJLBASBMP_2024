import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import argparse
import sys
import yaml

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

def downsample_single_seed(instance, seed, frac, downsampled_fastq_directory_specific):
    # for file in glob.glob(f"{instance.unzipped_fastq_directory}/*"):
    for file in glob.glob(f"{instance.original_fastq_directory}/*"):
        # cmd = ["seqtk", "sample", "-s", str(seed), file, str(frac)]
        output_file_path = os.path.join(downsampled_fastq_directory_specific, f"{os.path.basename(file)}")
        cmd = f"seqtk sample -s {seed} {file} {frac} | gzip > {output_file_path}"
        
        result = subprocess.run(cmd, stderr=subprocess.PIPE, shell=True, check=True)

        if result.returncode != 0:
            print(f"Error in processing {os.path.basename(file)} with seed {seed} and frac {frac}: {result.stderr.decode('utf-8')}")
        
        print(f"Finished downsampling file {os.path.basename(file)} with seed {seed} and frac {frac}")

def process_seed_frac_pair(instance, seed, frac):
    frac_str = str(frac).replace('.', '_')
    downsampled_fastq_directory_specific = os.path.join(
        instance.downsampled_fastq_directory,
        f"frac{frac_str}_seed{seed}")
    if not os.path.exists(downsampled_fastq_directory_specific):
        os.makedirs(downsampled_fastq_directory_specific)
    downsample_single_seed(instance, seed, frac, downsampled_fastq_directory_specific)
        
        
def downsample_fastqs_function(instance):
    os.chdir(instance.data_directory)
    for frac in instance.frac_list:
        for seed in instance.seed_list:
            if int(frac) < 1:
                print(f"downsampling with frac {frac}")
                process_seed_frac_pair(instance, seed, frac)
        print(f"Finished downsampling with frac {frac}")
        
    print("Finished downsampling!")
        
if __name__ == "__main__":
    from fastq_processor import FastqProcessor
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--yaml_path', default=f'{parent_path}/config.yaml', help='Path to the YAML configuration file.')
    args = parser.parse_args()
    
    with open(args.yaml_path, 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)
    
    print("entering!")

    fastq_processor.downsample_fastqs()
