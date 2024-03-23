import os
import sys
import yaml

import subprocess
import argparse
from datetime import datetime
import re
now = datetime.now()
date_directory_name = now.strftime('%y%m')

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

from scripts.organize_outputs import organize_output
from scripts.install_cellranger import install_cellranger_function

def ensure_correct_cellranger_version(instance):
    path_output_raw = subprocess.run("echo $PATH", shell=True, check=True, stdout=subprocess.PIPE, text=True)
    path_output = path_output_raw.stdout.strip()
    path_elements = path_output.split(':')
    new_path_elements = [element for element in path_elements if 'cellranger' not in element]

    desired_cellranger_path = f'{instance.package_path}/cellranger-{instance.cellranger_version}'
    
    # Define the regex pattern for the cellranger version
    current_cellranger_version_pattern = re.compile(r'cellranger-(\d+\.\d+\.\d+)')

    # Initialize variable to store the version
    current_cellranger_version = None
    all_cellranger_versions = []

    # Iterate through path elements to find the first cellranger version
    for element in path_elements:
        match = current_cellranger_version_pattern.search(element)
        if match:
            # Extract the version number
            if current_cellranger_version is None:
                current_cellranger_version = match.group(1)
            all_cellranger_versions.append(match.group(1))

    # Append the desired cellranger version at the beginning
    if current_cellranger_version != instance.cellranger_version:
        if instance.cellranger_version not in all_cellranger_versions:
            user_response = input(f"Cellranger {instance.cellranger_version} is not installed. Would you like to install it? (y/n): ")
            if user_response.lower() == 'y':
                install_cellranger_function(instance)
            else:
                print("Cellranger version in config.yaml not installed. To continue, please change config.yaml version or install this package")
                sys.exit("Installation aborted.")
        print(f"Successfully added cellranger version {instance.cellranger_version} as the current version in PATH.")
    else:
        print(f"Cellranger version {instance.cellranger_version} is correct")

    # Join the elements back into a PATH string
    new_path_elements.insert(0, desired_cellranger_path)
    updated_path = ':'.join(new_path_elements)
    os.environ["PATH"] = updated_path

def cellranger_install_transcriptome_function(instance):
    cellranger_reference_path = os.path.join(instance.reference_path, "cellranger", instance.cellranger_reference_folder_name)
    if not os.path.exists(cellranger_reference_path):
        os.makedirs(cellranger_reference_path)
    os.chdir(cellranger_reference_path)

    if not os.listdir(cellranger_reference_path):
        subprocess.run(f"wget {instance.cellranger_transcriptome_link}", shell=True, executable="/bin/bash")
        subprocess.run(f"tar -xzvf {instance.cellranger_transcriptome_link.split('/')[-1]}", shell=True, executable="/bin/bash")
        
def cellranger_count_function(instance, baseline):
    if baseline:
        seed_list = (0,)
        frac_list = (1.0,)
    else:
        seed_list = instance.seed_list
        frac_list = instance.frac_list

    if instance.cellranger_version != "":
        ensure_correct_cellranger_version(instance)

    cellranger_version_str = str(instance.cellranger_version).replace('.', '_')
    
    if instance.matrix_folder_name != "":
        out_folder_final = instance.matrix_folder_name
    else:
        out_folder_final = f"cellranger{cellranger_version_str}"

    cellranger_output_directory = os.path.join(instance.output_directory, out_folder_final)
    os.makedirs(cellranger_output_directory, exist_ok=True)
    os.chdir(cellranger_output_directory)
    cellranger_reference_path = os.path.join(instance.reference_path, "cellranger", instance.cellranger_reference_folder_name)

    transcriptome_path = None

    if os.path.exists(cellranger_reference_path):
        all_items = os.listdir(cellranger_reference_path)
        directories_only = [item for item in all_items if os.path.isdir(os.path.join(cellranger_reference_path, item))]
        for child_dir in directories_only:
            if child_dir.startswith("refdata"):
                transcriptome_path = os.path.join(cellranger_reference_path, child_dir)

    if transcriptome_path is None:
        raise Exception(f"Cellranger transcriptome not found. Please download it or check values in config.yaml")

    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            
            if baseline or int(frac) == 1:
                specific_fastq_directory = instance.original_fastq_directory
            else:
                specific_fastq_directory = os.path.join(instance.downsampled_fastq_directory, f"frac{frac_str}_seed{seed}")
            
            command = [
                "cellranger",
                "count",
                f"--fastqs={specific_fastq_directory}",
                f"--sample={instance.data_name}",
                "--localcores=8",
                "--localmem=64",
                f"--id=frac{frac_str}_seed{seed}",
                f"--transcriptome={transcriptome_path}"
            ]  # modify ID to change output folder name

            if instance.cellranger_include_introns == "true":
                if int(instance.cellranger_version[0]) < 7:
                    command.insert(2, "--include-introns")
            elif instance.cellranger_include_introns == "false":
                if int(instance.cellranger_version[0]) >= 7:
                    command.insert(2, "--include-introns=false")
             
            if instance.cellranger_expect_cells != "":        
                command.insert(2, f"--expect-cells={int(instance.cellranger_expect_cells)}")

            subprocess.run(command)

            print(f"cellranger count done for seed {seed} and frac {frac}")
            
            organize_output(instance.output_directory, seed, frac_str, matrix_source="cellranger", matrix_version = cellranger_version_str, matrix_folder_name = out_folder_final)

    print(f"cellranger count done!")


if __name__ == "__main__":    
    from fastq_processor import FastqProcessor
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', action='store_true', help='True if using baseline data (will override any seeds and counts accordingly), False if using downsampled fastq data')
    parser.add_argument('-c', '--download_transcriptome', action='store_true', help='Download reference transcriptome before running cellranger count')
    parser.add_argument('-y', '--yaml_path', default=f'{parent_path}/config.yaml', help='Path to the YAML configuration file.')
    args = parser.parse_args()
   
    with open(args.yaml_path, 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    if args.download_transcriptome:
        fastq_processor.cellranger_install_transcriptome()

    fastq_processor.cellranger_count(baseline=args.baseline)
