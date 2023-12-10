import subprocess
import os
from pygecko.parsers.utilities import list_files_and_directories





def msconvert(input_files: list|str, output_dir, format='mzML'):

    '''
    Takes in a list of input filenames and the output directory, converts the input files to the specified format and
    saves them in the output directory.

    Args:
        input_files (list|str): Input filename or list of input filenames.
        output_dir (str): Path to the output directory.
        format (str): Output format. Default: 'mzML'.
    '''

    # Define the path to the msconvert executable
    msconvert_path = r'C:\Users\felix\AppData\Local\Apps\ProteoWizard 3.0.23289.fd07aa9 64-bit\msconvert.exe'
    output_format = format

    # Check if the msconvert executable exists at the specified path
    if not os.path.exists(msconvert_path):
        print("msconvert executable not found at the specified path.")
        return

    # Change the working directory to the specified path
    #working_dir = r'C:\Users\felix\AppData\Local\Apps\ProteoWizard 3.0.23289.fd07aa9 64-bit'
    #os.chdir(working_dir)

    # Define the msconvert command
    msconvert_cmd = [
        msconvert_path,    # Path to the msconvert executable
        f'--{output_format}',         # Output format
        '-o', output_dir,  # Output directory
    ]

    if isinstance(input_files, str):
        input_files = list_files_and_directories(input_files)
    # Add the input filenames to the command
    msconvert_cmd.extend(input_files)

    try:
        # Execute the msconvert command
        subprocess.run(msconvert_cmd, check=True)
        print(f"Conversion to {output_format} format successful.")
    except subprocess.CalledProcessError as e:
        print(f"Error during conversion: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Provide a list of input filenames and the output directory
    input_files = ["C:/Users/felix/Documents/sciebo/AK/FKB-FA-005/FKB-FA-005_A.D", "C:Users/felix/Documents/sciebo/AK/FKB-FA-005/FKB-FA-005_B.D"]
    output_directory = "C:/Users/felix/Documents/sciebo/AK/Research/pyGECKO/pyGECKO/data"

    msconvert(input_files, output_directory)
    print()
