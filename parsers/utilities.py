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

if __name__ == "__main__":
    # Provide the path to the directory
    directory_path = r"C:\Users\felix\Documents\sciebo\AK\Research\pyGECKO\pyGECKO\data"

    # Call the function and store the result in a variable
    file_paths = list_files_and_directories(directory_path)

    # Print the result
    print(file_paths)
    print()