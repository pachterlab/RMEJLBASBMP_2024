import os
import sys
import yaml
import pkg_resources

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

import subprocess
import json
import shutil
from datetime import datetime
import argparse
from scripts.organize_outputs import organize_output
now = datetime.now()


def check_kb_version(target_version):
    version = pkg_resources.get_distribution("kb-python").version
    if version == target_version:
        print(f"kb version {version} matches config.yaml. Proceeding to run kb count")
    else:
        raise ValueError(f"kb target version = {target_version}, but actual kb version = {version}. Please check config.yaml and/or the python environment. kb-python can be installed with pip install kb-python gget ffq")

def custom_sort(filename):
    # Define order for file types
    file_type_order = {'R1': 0, 'R2': 1, 'I1': 2, 'I2': 3}

    # Split filename by '_' to extract file type and lane information
    parts = filename.split("_")

    # Extract lane number; assuming lane info is of the format 'L00X'
    lane = int(parts[-3][1:4])  # e.g., extracts '001' from 'L001'

    # Get the order value for the file type, e.g., 'R1'
    file_type = parts[-2].split(".")[0]  # e.g., extracts 'R1' from 'R1_001.fastq.gz'

    # Return a tuple: (file_type_order, lane_number)
    return (lane, file_type_order.get(file_type, 999))  # 999 as a default value for any unexpected file types



def kb_ref_function(instance):
    ensembl_path = os.path.join(instance.reference_path, "ensembl", str(instance.ensembl_release))
    
    if not os.path.exists(ensembl_path):
        os.makedirs(ensembl_path)
    os.chdir(ensembl_path)
    if not os.listdir(ensembl_path):
        gget_command = f"gget ref -o {ensembl_path}/info.json -r {instance.ensembl_release} -d -w dna,gtf homo_sapiens"
        print(gget_command)
        subprocess.run(gget_command, shell=True, executable="/bin/bash")
    
    kb_reference_path = os.path.join(instance.reference_path, "kb", instance.kb_reference_folder_name)

    if not os.path.exists(kb_reference_path):
        os.makedirs(kb_reference_path)
    os.chdir(kb_reference_path)
    if instance.kb_version != "":
        check_kb_version(instance.kb_version)
    
    if instance.kb_workflow == "nac":
        kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/fasta_spliced.fasta -f2 {kb_reference_path}/fasta_unspliced.fasta -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --workflow=nac {ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {ensembl_path}/Homo_sapiens.GRCh38.{instance.ensembl_release}.gtf.gz"
    else:   
        kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa --workflow=standard {ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {ensembl_path}/Homo_sapiens.GRCh38.{instance.ensembl_release}.gtf.gz"
    
    if instance.kb_d_list != "":
        position = "--workflow"
        position = kb_ref_command.find(position)
        kb_ref_command = kb_ref_command[:position] + f"--d-list {instance.kb_d_list} " + kb_ref_command[position:]

    if instance.kallisto_binary_path != "":
        position = f"{ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        position = kb_ref_command.find(position)
        kb_ref_command = kb_ref_command[:position] + f"--kallisto {instance.kallisto_binary_path} " + kb_ref_command[position:]

        result = subprocess.run([instance.kallisto_binary_path, "version"], capture_output=True, text=True)
        kallisto_version_output = result.stdout.strip()
        kallisto_major_version = kallisto_version_output.split('.')[1]
    else:
        kallisto_major_version = "50"

    if instance.bustools_binary_path != "":
        position = f"{ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        position = kb_ref_command.find(position)
        kb_ref_command = kb_ref_command[:position] + f"--bustools {instance.bustools_binary_path} " + kb_ref_command[position:]

    print(kb_ref_command)

    subprocess.run(kb_ref_command, shell=True, executable="/bin/bash")

    if int(kallisto_major_version) < 50:
        if instance.kallisto_binary_path != "":
            kallisto_binary_path = instance.kallisto_binary_path
        else:
            kallisto_binary_path_command = "kb info | grep 'kallisto:' | awk '{print $3}' | sed 's/[()]//g'"
            kallisto_binary_path = subprocess.run(kallisto_binary_path_command, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, text=True).stdout.strip()
        subprocess.run(f"{kallisto_binary_path} index -i {kb_reference_path}/index.idx -k 31 {kb_reference_path}/transcriptome.fa", shell=True, executable="/bin/bash")
    

        
def kb_count_function(instance, baseline, threads):
    if baseline:
        seed_list = (0,)
        frac_list = (1.0,)
    else:
        seed_list = instance.seed_list
        frac_list = instance.frac_list

    if instance.kb_version != "":
        check_kb_version(instance.kb_version)
    
    kb_version = pkg_resources.get_distribution("kb-python").version

    kb_version_major = int(kb_version.split('.')[1])

    kb_reference_path = os.path.join(instance.reference_path, "kb", instance.kb_reference_folder_name)

    if not os.path.exists(kb_reference_path) or (os.path.isdir(kb_reference_path) and not os.listdir(kb_reference_path)):
        raise Exception(f"{kb_reference_path} does not exist or is empty. Run kb ref or check config.yaml.")

    kb_version_str = str(kb_version).replace('.', '_')        
    
    if instance.matrix_folder_name != "":
        out_folder_final = instance.matrix_folder_name
    else:
        out_folder_final = f"kb{kb_version_str}"

    print("about to enter kb count")
    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            
            if baseline or int(frac) == 1:
                specific_fastq_directory = instance.original_fastq_directory
            else:
                specific_fastq_directory = os.path.join(instance.downsampled_fastq_directory, f"frac{frac_str}_seed{seed}")

            kb_version_str = str(kb_version).replace('.', '_')

            out_dir = os.path.join(instance.output_directory, out_folder_final, f"frac{frac_str}_seed{seed}")
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            # else:
            #     response = input(f"Warning: Output folder for frac{frac_str}_seed{seed} already exists in this project. Do you want to continue? (y/n) ")
            #     if response.lower() != "y":
            #         print("Operation aborted by the user.")
            #         sys.exit()
            
            os.chdir(specific_fastq_directory)

            print(f"Current directory: {specific_fastq_directory}")
            
            # Build the command
            kb_count_command = [
                'kb', 'count', 
                '-i', f'{kb_reference_path}/index.idx', 
                '-g', f'{kb_reference_path}/t2g.txt', 
                '-x', f'{instance.kb_sequencing_technology}', 
                '-o', f"{out_dir}",
                '-t', str(threads)
            ]

            if instance.kb_count_strandedness != "":
                kb_count_command.append(f"--strand {instance.kb_count_strandedness}")

            if instance.kb_workflow == "nac":
                new_arguments = ["--workflow=nac", "-c1", f"{kb_reference_path}/c_spliced.txt", "-c2", f"{kb_reference_path}/c_unspliced.txt", "--sum=total"]
                for argument in new_arguments:
                    kb_count_command.append(argument)

            else:
       	        kb_count_command.append("--workflow=standard")

            if instance.kallisto_binary_path != "":
                kb_count_command.append(f"--kallisto {instance.kallisto_binary_path}")
            
            if instance.bustools_binary_path != "":
                kb_count_command.append(f"--bustools {instance.bustools_binary_path}")

            # Add the fastq files
            fastq_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
            fastq_files = [filename for filename in os.listdir(specific_fastq_directory) if filename.endswith(fastq_extensions)]
            
            sorted_files = sorted(fastq_files, key=custom_sort)

            sorted_files = [file for file in sorted_files if "_I1_" not in file and "_I2_" not in file]

            # Append sorted filenames to kb_count_command
            for file in sorted_files:
                kb_count_command.append(file)
            
            print(' '.join(kb_count_command))

            kb_count_command = ' '.join(kb_count_command)
            
            # Run the command
            subprocess.run(kb_count_command, shell=True, executable="/bin/bash")

            # organize_output(instance.output_directory, seed, frac_str, matrix_source = "kb", matrix_version = kb_version_str, matrix_folder_name = out_folder_final)
            
            print(f"kb count finished for seed {seed} and frac {frac_str}")
    
    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            organize_output(instance.output_directory, seed, frac_str, matrix_source = "kb", matrix_version = kb_version_str, matrix_folder_name = out_folder_final)

    print("kb count finished")




    
            
if __name__ == "__main__":
    from fastq_processor import FastqProcessor
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', action='store_true', help='True if using baseline data (will override any seeds and counts accordingly), False if using downsampled fastq data')
    parser.add_argument('-t', '--threads', default='8', help='Number of threads for kb count')
    parser.add_argument('-r', '--run_ref', action='store_true', help='Run kb ref in addition to kb count')
    parser.add_argument('-y', '--yaml_path', default=f'{parent_path}/config.yaml', help='Path to the YAML configuration file.')
    args = parser.parse_args()
   
    with open(args.yaml_path, 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    if args.run_ref:
        fastq_processor.kb_ref()
    
    fastq_processor.kb_count(baseline=args.baseline, threads=args.threads)
