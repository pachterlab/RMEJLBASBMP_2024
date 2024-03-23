import os
import sys
import shutil
import subprocess
import yaml

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

def install_seqtk_function(instance):
    if not os.path.exists(os.path.join(instance.package_path, "seqtk")):
        subprocess.run(f"""
            cd {instance.package_path} &&
            git clone https://github.com/lh3/seqtk.git &&
            cd seqtk &&
            make &&
            chmod +x {instance.package_path}/seqtk/seqtk
            export PATH=$PATH:{instance.package_path}/seqtk
        """, shell=True, executable="/bin/bash")
    
    os.environ["PATH"] = f":{instance.package_path}/seqtk" + os.environ["PATH"]

    bashrc_content = f"""
    if [[ ":$PATH:" != *":{instance.package_path}/seqtk:"* ]]; then
        export PATH="{instance.package_path}/seqtk:$PATH"
    fi
    """
    
    bashrc_path = os.path.expanduser('~/.bashrc')

    with open(bashrc_path, "r") as file:
        existing_content = file.read()

    if bashrc_content.strip() not in existing_content:
        with open(f"{os.path.expanduser('~')}/.bashrc", "a") as f:
            f.write("\n" + bashrc_content)

def install_cellranger_function(instance):
    if not os.path.exists(os.path.join(instance.package_path, f"cellranger-{instance.cellranger_version}")):
        subprocess.run(f"""
            cd {instance.package_path} &&
            
            echo 'curl -o cellranger-{instance.cellranger_version}.tar.gz "{instance.cellranger_package_link}"' &&
            curl -o cellranger-{instance.cellranger_version}.tar.gz "{instance.cellranger_package_link}" &&
            
            echo 'tar -xzvf cellranger-{instance.cellranger_version}.tar.gz' &&
            tar -xzvf cellranger-{instance.cellranger_version}.tar.gz &&
            
            export PATH={instance.package_path}/cellranger-{instance.cellranger_version}:$PATH
        """, shell=True, executable="/bin/bash")
    os.environ["PATH"] = f":{instance.package_path}/cellranger-{instance.cellranger_version}" + os.environ["PATH"]
    
    bashrc_content = f"""
    if [[ ":$PATH:" != *":{instance.package_path}/cellranger-{instance.cellranger_version}:"* ]]; then
        export PATH="{instance.package_path}/cellranger-{instance.cellranger_version}:$PATH"
    fi
    """

    bashrc_path = os.path.expanduser('~/.bashrc')

    with open(bashrc_path, "r") as file:
        existing_content = file.read()

    if bashrc_content.strip() not in existing_content:
        with open(bashrc_path, "a") as f:
            f.write("\n" + bashrc_content)

if __name__ == "__main__":
    from fastq_processor import FastqProcessor
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)
    
    fastq_processor.install_cellranger()
