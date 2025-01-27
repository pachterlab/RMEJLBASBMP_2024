import os
import subprocess
import shutil

### BEGINNING OF ARGUMENTS ###
matrix_generation_base_dir = "RMEJLBASBMP_2024/matrix_generation/count_matrices_additional_datasets"
analysis_base_dir = "RMEJLBASBMP_2024/analysis"
remote = "user@x.x.x.x"
mouse_reference_folder = "reference/mouse"
human_reference_folder = "reference/human_kb28_kallisto50"
matrix_generation_local = False
kb_count_print_only = False
threads = 8
### END OF ARGUMENTS ###

sample_name_mouse = "neuron_10k_v3"
technology_mouse = "10xv3"
strand_mouse = None

sample_name_nsclc = "vdj_v1_hs_nsclc_multi_5gex_t_b"
technology_nsclc = "10xv1"
strand_nsclc = "unstranded"

aligner_name = "kb"
aligner_version = "0_28_0"
fraction = "1_0"
seed = "0"


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

def custom_sort(filepath):
    # Define order for file types
    file_type_order = {'R1': 0, 'R2': 1, 'I1': 2, 'I2': 3}

    # Split the filepath into parts by '/'
    path_parts = filepath.split("/")
    
    # Extract the parent folder (2nd to last part)
    parent_folder = path_parts[-2]

    # Extract the filename (last part of the path)
    filename = path_parts[-1]

    # Split filename by '_' to extract file type and lane information
    parts = filename.split("_")

    # Extract lane number; assuming lane info is of the format 'L00X'
    lane = int(parts[-3][1:4])  # e.g., extracts '001' from 'L001'

    # Get the order value for the file type, e.g., 'R1'
    file_type = parts[-2].split(".")[0]  # e.g., extracts 'R1' from 'R1_001.fastq.gz'

    # Return a tuple to sort by:
    # 1. Alphabetically by parent folder
    # 2. Numerically by lane
    # 3. Order of file type (R1, R2)
    return (parent_folder, lane, file_type_order.get(file_type, 999))

def run_kb_count(fastq_folder, reference_folder, kb_count_out = ".", technology = "10xv3", strand = None):
    fastq_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')

    # Find all FASTQ files in the folder and its subdirectories
    fastq_files = [
        os.path.join(root, filename)
        for root, _, files in os.walk(fastq_folder)
        for filename in files
        if filename.endswith(fastq_extensions)
    ]

    # remove the I1 and I2 files (don't pass into kb count)
    if technology != "10xv1":
        filtered_files = [
            f for f in fastq_files
            if not any(x in os.path.basename(f) for x in ['I1', 'I2'])
        ]
    else:
        filtered_files = fastq_files

    # sort by R1, R2
    sorted_files = sorted(filtered_files, key=custom_sort)
    
    # Build the command
    kb_count_command = [
        'kb', 'count', 
        '-i', f'{reference_folder}/index.idx', 
        '-g', f'{reference_folder}/t2g.txt', 
        '-x', technology, 
        '-o', kb_count_out,
        '-t', str(threads)
    ] + sorted_files

    # Insert ["strand", strand] after the technology element
    if strand is not None:
        kb_count_command = (
            kb_count_command[:7 + 1] +  # Up to and including 'technology' (the 7th element)
            ["strand", strand] +                      # Elements to insert
            kb_count_command[7 + 1:]   # Remaining elements
        )

    kb_count_command = ' '.join(kb_count_command)

    print(kb_count_command)

    if not kb_count_print_only:
        subprocess.run(kb_count_command, shell=True, executable="/bin/bash")



# Defining directories
data_dir = os.path.join(matrix_generation_base_dir, "data")
os.makedirs(data_dir, exist_ok=True)







# Mouse
mouse_data_dir = os.path.join(data_dir, f"{sample_name_mouse}_fastqs")
kb_count_out_mouse = os.path.join(mouse_data_dir, "kb_count_out")
src_matrix_path_mouse = os.path.join(kb_count_out_mouse, "counts_unfiltered")
analysis_folder_mouse = f"{analysis_base_dir}/count_matrix_collection/{sample_name_mouse}/{aligner_name}{aligner_version}/frac{fraction}_seed{seed}"

os.chdir(data_dir)
if not os.path.exists(mouse_data_dir) or not os.listdir(mouse_data_dir):
    download_mouse_command = "curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_fastqs.tar"
    extract_mouse_command = "tar -xvf neuron_10k_v3_fastqs.tar"
    subprocess.run(download_mouse_command, shell=True, executable="/bin/bash")
    subprocess.run(extract_mouse_command, shell=True, executable="/bin/bash")

mouse_fastq_folder = mouse_data_dir   # os.path.join(mouse_data_dir, "gex")

# make reference if needed
if not os.path.exists(mouse_reference_folder) or not os.listdir(mouse_reference_folder):
    os.makedirs(mouse_reference_folder, exist_ok=True)
    mouse_reference_download_command = f"kb ref -d mouse -i {mouse_reference_folder}/index.idx -g {mouse_reference_folder}/t2g.txt"
    subprocess.run(mouse_reference_download_command, shell=True, executable="/bin/bash")

# run kb count
os.makedirs(kb_count_out_mouse, exist_ok=True)
if not os.listdir(kb_count_out_mouse):
    run_kb_count(mouse_fastq_folder, mouse_reference_folder, kb_count_out_mouse, technology = technology_mouse, strand = strand_mouse)

if matrix_generation_local:
    copy_matrix_contents(src_matrix_path_mouse, analysis_folder_mouse)
else:
    mkdir_command_mouse = f"mkdir -p {analysis_folder_mouse}"
    download_command_mouse = f"scp -r {remote}:'{src_matrix_path_mouse}/*' {analysis_folder_mouse}/"

    print("Run the following commands locally:", mkdir_command_mouse, download_command_mouse, sep="\n")




# nsclc
nsclc_data_dir = os.path.join(data_dir, f"{sample_name_nsclc}_fastqs")
kb_count_out_nsclc = os.path.join(nsclc_data_dir, "kb_count_out")
src_matrix_path_nsclc = os.path.join(kb_count_out_nsclc, "counts_unfiltered")
analysis_folder_nsclc = f"{analysis_base_dir}/count_matrix_collection/{sample_name_nsclc}/{aligner_name}{aligner_version}/frac{fraction}_seed{seed}"

os.chdir(data_dir)
if not os.path.exists(nsclc_data_dir) or not os.listdir(nsclc_data_dir):
    download_nsclc_command = "curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/5.0.0/vdj_v1_hs_nsclc_multi_5gex_t_b/vdj_v1_hs_nsclc_multi_5gex_t_b_fastqs.tar"
    extract_nsclc_command = "tar -xvf vdj_v1_hs_nsclc_multi_5gex_t_b_fastqs.tar"
    # change directory name from "fastqs" to "vdj_v1_hs_nsclc_5gex_fastqs"
    shutil.move("fastqs", "vdj_v1_hs_nsclc_multi_5gex_t_b_fastqs")
    subprocess.run(download_nsclc_command, shell=True, executable="/bin/bash")
    subprocess.run(extract_nsclc_command, shell=True, executable="/bin/bash")
    

nsclc_fastq_folder = os.path.join(nsclc_data_dir, "vdj_v1_hs_nsclc_5gex_fastqs")  # nsclc_data_dir

# make reference if needed
if not os.path.exists(human_reference_folder) or not os.listdir(human_reference_folder):
    os.makedirs(human_reference_folder, exist_ok=True)
    nsclc_reference_download_command = f"kb ref -d human -i {human_reference_folder}/index.idx -g {human_reference_folder}/t2g.txt"
    subprocess.run(nsclc_reference_download_command, shell=True, executable="/bin/bash")

# run kb count
os.makedirs(kb_count_out_nsclc, exist_ok=True)
if not os.listdir(kb_count_out_nsclc):
    run_kb_count(nsclc_fastq_folder, human_reference_folder, kb_count_out_nsclc, technology = technology_nsclc, strand = strand_nsclc)


if matrix_generation_local:
    copy_matrix_contents(src_matrix_path_nsclc, analysis_folder_nsclc)
else:    
    mkdir_command_nsclc = f"mkdir -p {analysis_folder_nsclc}"
    download_command_nsclc = f"scp -r {remote}:'{src_matrix_path_nsclc}/*' {analysis_folder_nsclc}/"

if not matrix_generation_local:
    print("Run the following commands locally:", mkdir_command_mouse, download_command_mouse, mkdir_command_nsclc, download_command_nsclc, sep="\n")