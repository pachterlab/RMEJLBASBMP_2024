import os
import shutil
import glob
import subprocess
import argparse
import yaml
import sys

current_path = os.path.dirname(os.path.abspath(__file__))  # gets the current directory of the file that this script is in
sys.path.append(current_path) if current_path not in sys.path else None

from fastq_processor import FastqProcessor 


def arg_parser_setup():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-y', '--yaml_path', default='config.yaml', help='Path to the YAML configuration file.')
    
    parser.add_argument(
        "-p",
        "--project_setup",
        action="store_true",
        help="Enable project setup - install seqtk, install cellranger command line, organizes folder structure",
    )
    parser.add_argument(
        "-c",
        "--download_cellranger_transcriptome",
        action="store_true",
        help="download cellranger transcriptome",
    )
    parser.add_argument(
        "-k", "--download_kb_index", action="store_true", help="download kb index"
    )
    parser.add_argument(
        "-q", "--download_fastqs", action="store_true", help="download fastqs"
    )
    parser.add_argument(
        "-f",
        "--download_filtered_feature_bc_matrix",
        action="store_true",
        help="download filtered feature matrix",
    )
    
    parser.add_argument("-a", "--all", action="store_true", help="Enable all options (good for starting from scratch)")

    args = parser.parse_args()

    if args.all:
        args.project_setup = True
        args.download_cellranger_transcriptome = True
        args.download_kb_index = True
        args.download_fastqs = True
        args.download_filtered_feature_bc_matrix = True

    return args


def main(args):
    yaml_file_path = args.yaml_path
    
    with open(yaml_file_path, 'r') as file:
        config = yaml.safe_load(file)
    
    fastq_processor = FastqProcessor(**config)

    if args.project_setup:
        fastq_processor.project_setup()

    if args.download_filtered_feature_bc_matrix:
        os.chdir(os.path.join(fastq_processor.output_directory, "10x_genomics_matrices"))
        subprocess.run(f"wget {fastq_processor.filtered_feature_bc_matrix_link}", shell=True, executable="/bin/bash")
    
    if args.download_fastqs and fastq_processor.fastq_link is not None:
        os.chdir(fastq_processor.data_directory)
        
        if os.path.exists(fastq_processor.original_fastq_directory):
            if not os.listdir(fastq_processor.original_fastq_directory):
                subprocess.run(f"wget {fastq_processor.fastq_link}", shell=True, executable="/bin/bash")
                subprocess.run(f"tar -xvf {fastq_processor.data_name}_fastqs.tar", shell=True, executable="/bin/bash")
        else:
            subprocess.run(f"wget {fastq_processor.fastq_link}", shell=True, executable="/bin/bash")
            subprocess.run(f"tar -xvf {fastq_processor.data_name}_fastqs.tar", shell=True, executable="/bin/bash")
    
    if fastq_processor.frac_list:
        fastq_processor.downsample_fastqs()
    # fastq_processor.read_counts(fastq_processor.data_name, fastq_processor.data_directory, fastq_processor.output_directory, fastq_processor.seed_list, fastq_processor.frac_list)

    if fastq_processor.count_matrix_generation_method == "kb" or fastq_processor.count_matrix_generation_method == "both":
        if args.download_kb_index:
            fastq_processor.kb_ref()

        # kb count on full size data
        fastq_processor.kb_count(baseline=True)

        # kb count on downsampled data
        if fastq_processor.frac_list:
            fastq_processor.kb_count(baseline=False)

    if fastq_processor.count_matrix_generation_method == "cellranger" or fastq_processor.count_matrix_generation_method == "both":
        if args.download_cellranger_transcriptome:
            fastq_processor.cellranger_install_transcriptome()
        fastq_processor.cellranger_count(baseline=True)
        if fastq_processor.frac_list:
            fastq_processor.cellranger_count(baseline=False)

    print("Main ran successfully. Please find final data is in count_matrix_collection/<project_name>/<matrix_source>_filtered_feature_bc_matrix_collection")


if __name__ == "__main__":
    args = arg_parser_setup()
    
    main(args)
