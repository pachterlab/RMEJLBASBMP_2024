# to just make conda environment
conda env create -f environment.yml

# select appropriate config.yaml file and change parameters as desired

# minimal setup (requires fastq files, reference transcriptome (cellranger) or index/t2g (kb) files in reference_files/, and cellranger/kb-python and seqtk installed, all in the structure described below):
python3 main.py

# include all data downloads (warning: will download large fastq files - see more below):
python3 main.py -a

# can also pick and choose which functions to run by running the appropriate file in src


# To run analysis:
docker run -d -p 8787:8787 -e PASSWORD=yp -v /path/to/project:/workspace josephrich98/seurat_vs_scanpy:1.0
# Then go to localhost:8787 in your browser and login with username rstudio and password yp


### note: for baseline (no downsampling), notation is seed=0, frac=1


### main.py without any flags does not download any large files, but only works with the following file structure. Please either manually ensure the following file structure (especially noting the data folder and folders outside project_root) or run main.py with all flags included:

~/<root>/  
├── opt/  
│   ├── cellranger-X.X.X/  
│   └── seqtk  
├── reference_files/  
│   ├── cellranger/  
│   │   └── <date_of_install>/  
│   │       └── grch38_transcriptome/  
│   │           └── refdata-gex-GRCh38-2020-A/  
│   │               └── etc.  
│   └── kb/  
│       └── <date_of_install>/
│           ├── index.idx/
│           ├── t2g.txt
│           └── transcriptome.fa
│
└── matrix_generation/
    ├── README.md
    ├── main.py
    ├── conda_environments
    ├── scripts/
    ├── src/
    ├── tests/
    ├── data/
    │   └── <data_name>/
    │       ├── fastqs/
    │       │   ├── file1.fastq.gz
    │       │   ├── file2.fastq.gz
    │       │   └── etc.
    │       └── downsampled_fastqs/
    │           ├── file1.fastq
    │           ├── file2.fastq
    │           └── etc.
    └── count_matrix_collection/
        └── <data_name>/
            ├── read_counts/
            └── <matrix_source>/
                ├── filtered_feature_bc_matrix_collection
                │   └──<FINAL_MATRICES_HERE>
                └── <run_name>/
                    └── <additional_parent_file(s)>
                        ├── <matrix>.mtx
                        ├── <genes>.tsv
                        └── <barcodes>.tsv



Commands listed if preferred to run manually:
- installs (seqtk, cellranger, kb-python)
cd /home/jrich/Desktop/opt && git clone https://github.com/lh3/seqtk.git && cd seqtk && make && chmod +x /home/jrich/Desktop/opt/seqtk/seqtk export PATH=$PATH:/home/jrich/Desktop/opt/seqtk

cd /home/jrich/Desktop/opt &&
curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1695707214&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTU3MDcyMTR9fX1dfQ__&Signature=ArH5uCG~D3PL3JkIUrhZXbupvWHKOzcJ0oFW-cIQinJBFlTgZNOA-9zBZQ8vQWw5fbeoQAdtLySEcpAn8Cb0lPDheoXcCexePGHzrRg65bEpv1-A6tB8pJ9Tbyr8l4G4CKtn1zl5ov3X5PGY2~eM1iMplh3EfK4ES6t9Zpju~uxnm0Dybj0m59eZMl15L-nbb1Q2Vy4jKZYb5BMCtNDulVvbabVFbx80UJFbS~W-LsCk3okWKoZ3ESal-TOcRJVvYgbYtPp~Lgii9ZuizSdjla3aL7Tg3LwIfNFvGfzBKywNfrVv~7HkLnjkeDNgOb1x54jgWYl0Gski6S~cmEpfGA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && tar -xzvf cellranger-7.1.0.tar.gz && export PATH=$PATH:{package_path}/cellranger

conda create python=3.11
KB28 (kallisto v50): pip install kb-python==0.28.0 gget==0.27.9 pyyaml==6.0.1
KB24 (kallisto v46): pip install kb-python==0.24.4 gget==0.27.9 pyyaml==6.0.1 ngs-tools==1.8.5


- downloads (reference transcriptome/index files, fastq files)
cd /home/jrich/Desktop/reference_files/cellranger/<date_directory_name> && wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz && tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz && rm refdata-gex-GRCh38-2020-A.tar.gz

cd /home/jrich/Desktop/downsample_setup/data/<project_name> && wget <fastq_link> && tar -xvf <data_name>_fastqs.tar

- downsample fastq files (seqtk)
seqtk sample -s<seed> <fastqfile> <fastqfile> <frac> | gzip > {output_file_path}"

- compute downsampled read counts
python3 src/downsampled_read_counts.py <project_name> <frac>

- IF KB run kb-ref
kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa $(gget ref --ftp -w dna,gtf homo_sapiens)

- run cellranger/kb count 
cd ~/Desktop/data/<project_name>
kb count -i index.idx -g t2g.txt -x 10xv3 -o ~/Desktop/output/<project_name>/kb FASTQ1 FASTQ2 ETC
cellranger count --id=seed{seed}_frac{int(frac*100)} --transcriptome={transcriptome_path} --fastqs={fastq_directory} --sample={instance.data_name} --localcores=8 --localmem=64

