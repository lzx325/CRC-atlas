import os
import requests
import hashlib
from pathlib import Path
from tqdm import tqdm

def download_if_not_exist(fp, url, md5sum=None):
    """
    Download file from the given URL to the given fp if it doesn't exist.

    :param url: The URL of the file to download.
    :param fp: The local fp to save the file as.
    """

    if not os.path.exists(fp) or (md5sum is not None and compute_md5(fp) != md5sum):
        f_dir=Path(fp).parent
        os.makedirs(f_dir,exist_ok=True)

        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024  # 1 Kbyte
        t = tqdm(total=total_size, unit='B', unit_scale=True, desc=fp)

        with open(fp, 'wb') as f:
            for data in response.iter_content(block_size):
                t.update(len(data))
                f.write(data)

        t.close()
        if total_size != 0 and t.n != total_size:
            raise ValueError("ERROR, something went wrong")
        if md5sum is not None and compute_md5(fp) != md5sum:
            raise ValueError(f"MD5 checksum mismatch for file {fp}!")
        
        print(f"{fp} downloaded successfully!")

    else:
        if md5sum is not None:
            print(f"{fp} already exists and passed md5 check.")
        else:
            print(f"{fp} already exists.")

def compute_md5(filename):
    """
    Compute the MD5 hash of a given file.

    :param filename: Path to the file.
    :return: Hexadecimal MD5 hash string.
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as f:
        # Read and update hash in chunks to save memory
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()