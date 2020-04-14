import os
import pathlib
import time
import urllib.request

def download(url):
    """Download and return the content of a URL."""
    print(f"Downloading {url}... ", end="", flush=True)
    req = urllib.request.urlopen(url)
    data = req.read()
    print("Done.", flush=True)
    return data


def download_and_save(url, path, cache_duration=1000000000, load=True):
    """Download the URL, store to a file, and return its content.

    Arguments:
        url: URL to download.
        path: Target file path.
        cache_duration: (optional) Reload if the file on disk is older than the given duration in seconds.
        load: Should the file be loaded in memory? If not, `None` is returned.
    """
    path = pathlib.Path(path)
    try:
        if time.time() - path.lstat().st_mtime <= cache_duration:
            if load:
                with open(path, 'rb') as f:
                    return f.read()
            elif os.path.exists(path):
                return None
    except FileNotFoundError:
        pass

    data = download(url)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'wb') as f:
        f.write(data)
    if load:
        return data
    else:
        return None
