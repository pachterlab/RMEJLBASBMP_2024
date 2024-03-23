import os 
import glob
import argparse
import sys
import yaml
import gzip

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

def read_counts_function(instance):
    read_counts = {}
    read_output_directory = os.path.join(instance.output_directory, "read_counts")
    if not os.path.exists(read_output_directory):
        os.makedirs(read_output_directory)
    read_output_file = os.path.join(read_output_directory, "baseline.txt")
    if not os.path.exists(read_output_file):    
        # file_list = glob.glob("/home/jrich/Desktop/downsample_setup/tests/*")
        file_list = glob.glob(f"{instance.original_fastq_directory}/*")
        
        file_list = [file for file in file_list if file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
        
        for file in file_list:
            base_name = os.path.basename(file)
            count = 0
            
            if file.endswith('.gz'):
                with gzip.open(file, 'rt') as f:  # 'rt' mode for reading as text
                    for _ in f:
                        count += 1

            else:
                with open(file, 'r') as f:
                    for _ in f:
                        count += 1

            read_counts[base_name] = int(count/4)
        
        with open(f"{read_output_directory}/baseline.txt", "w") as f:
            for k, v in read_counts.items():
                f.write(f"{k}: {v}\n")
                
    else:
        # # If reading counts from file
        # Open the file for reading
        with open(f"{read_output_file}", "r") as f:
            for line in f:
                # Remove any leading and trailing whitespaces (like newline characters)
                line = line.strip()

                # Split the line by ': ' into key and value parts
                key, value = line.split(': ')

                # Add the key-value pair to the dictionary
                read_counts[key] = int(value)  # Convert the value to integer before storing

    print(read_counts)  # Print the dictionary to verify it's as expected


    for frac in instance.frac_list:
        for seed in instance.seed_list:
            frac_str = str(frac).replace('.', '_')
            new_read_counts = {}
            downsampled_specific_path = os.path.join(instance.downsampled_fastq_directory, f"frac{frac_str}_seed{seed}")

            ## if I want to test a whole directory
            file_list = glob.glob(f"{downsampled_specific_path}/*")

            file_list = [file for file in file_list if file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
            
            for file in file_list:
                base_name = os.path.basename(file)
                count = 0

                if file.endswith('.gz'):
                    with gzip.open(file, 'rt') as f:  # 'rt' mode for reading as text
                        for _ in f:
                            count += 1

                else:
                    with open(file, 'r') as f:
                        for _ in f:
                            count += 1

                print(f"SEED {seed}, FRAC {frac}, FILE {base_name}: ", int(count/4))
                new_read_counts[base_name] = int(count/4)
                if not ((frac - 0.05) < new_read_counts[base_name] / read_counts[os.path.basename(file)] < (frac + 0.05)):
                    print(f"SEED {seed}, FRAC {frac}, FILE {base_name}: IS BAD - should be near {read_counts[os.path.basename(file)]*frac}, file is {new_read_counts[base_name]}")
            
            with open(f"{read_output_directory}/seed{seed}_frac{frac_str}.txt", "w") as f:
                for k, v in new_read_counts.items():
                    f.write(f"{k}: {v}\n")
                    if not ((frac - 0.05) < v / read_counts[k] < (frac + 0.05)):
                        f.write(f"SEED {seed}, FRAC {frac}, FILE {k}: IS BAD - should be near {read_counts[k]*frac}, file is {v}\n")

if __name__ == "__main__":
    from fastq_processor import FastqProcessor
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--yaml_path', default=f'{parent_path}/config.yaml', help='Path to the YAML configuration file.')
    args = parser.parse_args()
    
    with open(args.yaml_path, 'r') as file:
        config = yaml.safe_load(file)
    
    fastq_processor = FastqProcessor(**config)
    
    fastq_processor.read_counts()
