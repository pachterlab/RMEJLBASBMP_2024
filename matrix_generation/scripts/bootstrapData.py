import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import random
import os
import pyfastx
from tqdm import tqdm

### BEGINNING OF ARGUMENTS ###
file_list = [("SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L002_R1_001.fastq.gz",
             "SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L002_R2_001.fastq.gz"),
             ("SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L003_R1_001.fastq.gz",
             "SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L003_R2_001.fastq.gz"),
             ("SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L004_R1_001.fastq.gz",
             "SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L004_R2_001.fastq.gz")]  # either a list of strings (single-end) or a list of 2-tuples of strings (paired-end)

seed_list = [42, 43, 44, 45, 46]
threads = 1
bootstrapped_output_directory = "bootstrapped_fastqs_pbmc10k"
gzip_output = True
use_buffer = True
batch_size = 200000  # Number of reads to write in each batch; only if use_buffer is True

fastq_to_length_dict = {'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L002_R1_001.fastq.gz': 631503136, 'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L002_R2_001.fastq.gz': 631503136, 'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L003_R1_001.fastq.gz': 632331834, 'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L003_R2_001.fastq.gz': 632331834, 'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L004_R1_001.fastq.gz': 625535751, 'SC3_v3_NextGem_SI_PBMC_10K_fastqs/SC3_v3_NextGem_SI_PBMC_10K_S1_L005_R2_001.fastq.gz': 625535751}
### END OF ARGUMENTS ###

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, open_func, write_mode, seed = None):
    input_fastq_file_name = str(input_fastq).split()[-1]
    buffer = []  # Temporary storage for the batch
    if use_buffer:
        with open_func(output_path, write_mode) as f:
            for i, (name, seq, qual) in enumerate(tqdm(input_fastq, desc=f"Iterating through seed {seed}, file {input_fastq_file_name}", unit="read", total=total_reads)):
                # Add the FASTQ entry to the buffer
                buffer.extend([f"@{name}\n{seq}\n+\n{qual}\n"] * occurrence_list[i])
                
                # If the buffer reaches the batch size, write all at once and clear the buffer
                if (i + 1) % batch_size == 0:
                    f.writelines(buffer)
                    buffer.clear()  # Clear the buffer after writing
            
            # Write any remaining entries in the buffer
            if buffer:
                f.writelines(buffer)
                buffer.clear()
    else:
        with open_func(output_path, write_mode) as f:
            for i, (name, seq, qual) in enumerate(tqdm(input_fastq, desc=f"Iterating through seed {seed}, file {input_fastq_file_name}", unit="read", total=total_reads)):
                f.writelines(
                    f"@{name}\n{seq}\n+\n{qual}\n" * occurrence_list[i]
                )

def bootstrap_single_file(file = None, file1 = None, file2 = None, gzip_output = None, output_path = None, output_path1 = None, output_path2 = None, seed = None):
    # args = parse_arguments()

    if gzip_output:
        open_func = gzip.open
        write_mode = "wt"
    else:
        open_func = open
        write_mode = "w"

    # Create output directory if it doesn't exist
    output_path_args = [output_path, output_path1, output_path2]
    for output_path_arg in output_path_args:
        if output_path_arg and os.path.dirname(output_path_arg):
            os.makedirs(os.path.dirname(output_path_arg), exist_ok=True)

    if file:
        if not output_path:
            output_path = file.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path.endswith(".gz"):
            output_path += ".gz"

        # Single FASTQ file
        input_fastq_read_only = pyfastx.Fastx(file)
                
        print(f"Calculating total reads and determining random indices for seed {seed}, file {file}")
        total_reads = fastq_to_length_dict[file]
        random.seed(seed)
        random_indices = random.choices(range(total_reads), k=total_reads)

        # Initialize a list with zeros
        occurrence_list = [0] * total_reads

        # Count occurrences (I don't use a counter in order to save memory, as a counter is essentially a dictionary)
        for index in tqdm(random_indices, desc=f"Counting occurrences for seed {seed}, file {file}", unit="read", total=total_reads):
            occurrence_list[index] += 1

        del random_indices

        # write fastq
        write_fastq(input_fastq = input_fastq_read_only, output_path = output_path, occurrence_list = occurrence_list, total_reads = total_reads, open_func = open_func, write_mode = write_mode, seed = seed)

    elif file1 and file2:
        if not output_path1:
            output_path1 = file1.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path1.endswith(".gz"):
            output_path1 += ".gz"
        if not output_path2:
            output_path2 = file.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path2.endswith(".gz"):
            output_path2 += ".gz"

        input_fastq1_read_only = pyfastx.Fastx(file1)
        input_fastq2_read_only = pyfastx.Fastx(file2)
                
        # Paired-end FASTQ files
        print(f"Calculating total reads and determining random indices for seed {seed}, file {file1}")
        total_reads = fastq_to_length_dict[file1]
        random.seed(seed)
        random_indices = random.choices(range(total_reads), k=total_reads)

        # Initialize a list with zeros
        occurrence_list = [0] * total_reads

        # Count occurrences (I don't use a counter in order to save memory, as a counter is essentially a dictionary)
        for index in tqdm(random_indices, desc=f"Counting occurrences for seed {seed}, file {file1}", unit="read", total=total_reads):
            occurrence_list[index] += 1

        del random_indices

        # write fastqs
        write_fastq(input_fastq = input_fastq1_read_only, output_path = output_path1, occurrence_list = occurrence_list, total_reads = total_reads, open_func = open_func, write_mode = write_mode, seed = seed)
        write_fastq(input_fastq = input_fastq2_read_only, output_path = output_path2, occurrence_list = occurrence_list, total_reads = total_reads, open_func = open_func, write_mode = write_mode, seed = seed)

    else:
        raise ValueError("You must provide either a single FASTQ file or paired-end FASTQ files.")


def process_seed_and_file(seed, file):
    # print(f"seed {seed}, file {file}")
    if isinstance(file, tuple):
        assert len(file) == 2, "Paired-end FASTQ files must be a 2-tuple."
        output_path1 = os.path.join(bootstrapped_output_directory, f"frac1_0_seed{seed}", os.path.basename(file[0]))
        output_path2 = os.path.join(bootstrapped_output_directory, f"frac1_0_seed{seed}", os.path.basename(file[1]))
        bootstrap_single_file(file1 = file[0], file2 = file[1], gzip_output = gzip_output, output_path = None, output_path1 = output_path1, output_path2 = output_path2, seed = seed)
    elif isinstance(file, str):
        output_path = os.path.join(bootstrapped_output_directory, f"frac1_0_seed{seed}", os.path.basename(file))
        bootstrap_single_file(file = file, gzip_output = gzip_output, output_path = output_path, seed = seed)

def bootstrap_multiple_files(file_list, seed_list, threads=2):
    # Use ThreadPoolExecutor to process seeds in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit tasks for all combinations of seeds and files
        futures = [
            executor.submit(process_seed_and_file, seed, file)
            for seed in seed_list
            for file in file_list
        ]

def count_reads(filepath):
    fastq_file = pyfastx.Fastx(filepath)
    num_reads = sum(1 for _ in fastq_file)
    return num_reads

if 'fastq_to_length_dict' not in globals():
    fastq_to_length_dict = {}
def make_fastq_to_length_dict(file_list):
    global fastq_to_length_dict
    for file in file_list:
        if isinstance(file, tuple):
            assert len(file) == 2, "Paired-end FASTQ files must be a 2-tuple."
            if file[0] in fastq_to_length_dict and file[1] in fastq_to_length_dict:
                continue
            print(f"Counting {file[0]}")
            count = count_reads(file[0])
            fastq_to_length_dict[file[0]] = count
            fastq_to_length_dict[file[1]] = count
        elif isinstance(file, str):
            if file in fastq_to_length_dict:
                continue
            print(f"Counting {file}")
            count = count_reads(file)
            fastq_to_length_dict[file] = count
    print(fastq_to_length_dict)


if __name__ == "__main__":
    make_fastq_to_length_dict(file_list)
    bootstrap_multiple_files(file_list=file_list, seed_list=seed_list, threads=threads)

    