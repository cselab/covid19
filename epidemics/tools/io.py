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


def download_and_save(url, path, cache_duration=0):
    """Download the URL, store to a file, and return its content.

    Arguments:
        url: URL to download.
        path: Target file path.
        cache_duration: (optional) Reload if the file on disk is older than the given duration in seconds.
    """
    path = pathlib.Path(path)
    try:
        if time.time() - path.lstat().st_mtime <= cache_duration:
            with open(path, 'rb') as f:
                return f.read()
    except FileNotFoundError:
        pass

    data = download(url)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'wb') as f:
        f.write(data)
    return data
