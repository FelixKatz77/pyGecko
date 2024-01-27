import sys, os
from pathlib import Path


def list_files_and_directories(directory_path):
    # Create a Path object for the directory
    directory = Path(directory_path)

    # Initialize an empty list to store absolute paths
    file_paths = []

    # Use the iterdir() method to list files and directories in the specified directory
    for item in directory.iterdir():
        file_paths.append(str(item))

    return file_paths

def find_directories_with_extension(directory_path, extension):
    directory_path = Path(directory_path)
    if not directory_path.is_dir():
        print(f"{directory_path} is not a valid directory.")
        return []

    matching_directories = [entry for entry in directory_path.iterdir() if entry.is_dir() and entry.suffix == extension]

    return matching_directories

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
