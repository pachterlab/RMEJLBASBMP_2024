import os
import shutil
import subprocess

### BEGINNING OF ARGUMENTS ###
bootstrapped_output_directory = "bootstrapped_fastqs_pbmc10k"  # from bootstrapData.py
kb_count_output_parent_directory = "bootstrapped_count_matrices"
threads = 8
remote = "user@x.x.x.x"
kb_reference_path = "reference/human_kb28_kallisto50"
analysis_base_dir = "RMEJLBASBMP_2024/analysis"
kb_count_print_only = False
matrix_generation_local = False
### END OF ARGUMENTS ###

def copy_matrix_contents(src_matrix_path, dest_matrix_path):
    os.makedirs(dest_matrix_path, exist_ok=True)

    # Copy the contents of myfolder to mydestination
    for item in os.listdir(src_matrix_path):
        source = os.path.join(src_matrix_path, item)
        destination = os.path.join(dest_matrix_path, item)

        if os.path.isdir(source):
            shutil.copytree(source, destination, dirs_exist_ok=True)  # Copy directory
        else:
            shutil.copy2(source, destination)  # Copy file

# see download_and_align_additional_datasets for most updated custom_sort and related behavior (that generalized beyond this PBMC 10k dataset)
def custom_sort(filename):
    # Define order for file types
    file_type_order = {'R1': 0, 'R2': 1}

    # Split filename by '_' to extract file type and lane information
    parts = filename.split("_")

    # Extract lane number; assuming lane info is of the format 'L00X'
    lane = int(parts[-3][1:4])  # e.g., extracts '001' from 'L001'

    # Get the order value for the file type, e.g., 'R1'
    file_type = parts[-2].split(".")[0]  # e.g., extracts 'R1' from 'R1_001.fastq.gz'

    # Return a tuple: (lane, file_type_order), with 999 as a default for unexpected types
    return (lane, file_type_order.get(file_type, 999))

def run_kb_count_on_bootstrapped_data():
    fastq_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')

    # assumes directory structure is {bootstrapped_output_directory}/frac_1_0_seed{seed}/samplename_R1_L001_etc.fastq.gz
    for folder_name in os.listdir(bootstrapped_output_directory):
        folder_path = os.path.join(bootstrapped_output_directory, folder_name)
        if os.path.isdir(folder_path):  # Check if it's a directory
            fastq_files = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if filename.endswith(fastq_extensions)]
            filtered_files = [f for f in fastq_files if not any(x in f for x in ['I1', 'I2'])]  # remove the I1 and I2 files (don't pass into kb count)
            sorted_files = sorted(filtered_files, key=custom_sort)

            kb_count_out = os.path.join(kb_count_output_parent_directory, folder_name)

            # Build the command
            kb_count_command = [
                'kb', 'count', 
                '-i', f'{kb_reference_path}/index.idx', 
                '-g', f'{kb_reference_path}/t2g.txt', 
                '-x', '10xv3', 
                '-o', kb_count_out,
                '-t', str(threads)
            ] + sorted_files

            kb_count_command = ' '.join(kb_count_command)

            print(kb_count_command)

            if not kb_count_print_only:
                subprocess.run(kb_count_command, shell=True, executable="/bin/bash")

if __name__ == "__main__":
    run_kb_count_on_bootstrapped_data()

    analysis_folder_pbmc = f"{analysis_base_dir}/count_matrix_collection/SC3_v3_NextGem_SI_PBMC_10K/kb0_28_0"

    if matrix_generation_local:
        copy_matrix_contents(kb_count_output_parent_directory, analysis_folder_pbmc)
        print(f"Final output in {analysis_folder_pbmc}")
    else:    
        mkdir_command_pbmc = f"mkdir -p {analysis_folder_pbmc}"
        download_command_pbmc = f"scp -r {remote}:'{kb_count_output_parent_directory}/*' {analysis_folder_pbmc}/"
        print("Run the following commands locally:", mkdir_command_pbmc, download_command_pbmc, sep="\n")
    # then must move the downloaded files from parent_folder/counts_unfiltered/matrix_files to parent_folder/matrix_files

    