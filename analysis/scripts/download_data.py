import os
import tarfile
import requests
import json
from tqdm import tqdm

def download_file(doi, filename, data_path_root, overwrite = False):
    full_path = os.path.join(data_path_root, filename)
    if os.path.exists(full_path) and not overwrite:
        # raise FileExistsError(f"File {full_path} already exists.")
        print("File exists. Skipping download.")
        return full_path

    if doi.startswith('https://doi.org/'):
        doi = doi.replace('https://doi.org/', '')
    url = 'https://api.datacite.org/dois/'+doi+'/media'
    r = requests.get(url).json()
    found = False
    for item in r['data']:
        url = item['attributes']['url']
        # Extract the last part of the URL after the last slash
        url_filename = url.split('/')[-1]
        if url_filename == filename:
            netcdf_url = url
            found = True
            break  # Exit the loop since we've found the matching filename
    if not found:
        print("Error: No matching filename found in the data.")
        return None
    r = requests.get(netcdf_url,stream=True)
    
    os.makedirs(data_path_root, exist_ok=True)
    
    #Download file with progress bar
    if r.status_code == 403:
        print("File Unavailable")
        return None
    if 'content-length' not in r.headers:
        print("Did not get file")
        return None
    else:
        with open(full_path, 'wb') as f:
            total_length = int(r.headers.get('content-length'))
            for chunk in tqdm(r.iter_content(chunk_size=1024), total=total_length/1024, unit="KB"):
                if chunk:
                    f.write(chunk)
        return full_path

def get_true_folder_name_from_tarball(tar_path):
    with tarfile.open(tar_path, "r:gz") as tar:
        # Read first member
        first_member = tar.getmembers()[0]
        top = first_member.name.split("/")[0]
    return top

def download_and_extract(doi, filename, data_path_root, overwrite = False):
    full_path = download_file(doi, filename, data_path_root, overwrite=overwrite)
    if not os.path.exists(full_path):
        raise ValueError("Issue downloading")

    file_ext = filename.split(".", 1)[1]
    if file_ext == 'tar.gz':
        true_folder_name = get_true_folder_name_from_tarball(full_path)
        extracted_folder_name = full_path.replace(filename, true_folder_name)
        if not os.path.exists(extracted_folder_name):
            try:
                with tarfile.open(full_path, 'r:gz') as tar:
                    tar.extractall(path=full_path)
                print(f"Extraction complete: {filename}. Contents are in {full_path}")
            except Exception as e:
                print(f"Error extracting the file: {e}")
    return full_path
