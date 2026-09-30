import logging
import subprocess
import os
import shutil
from pygecko.parsers.utilities import list_files_and_directories

logger = logging.getLogger(__name__)


def find_msconvert() -> str | None:

    '''
    Returns the path to the msConvert executable: the PYGECKO_MSCONVERT environment variable if set, otherwise
    an `msconvert` found on PATH. None if neither points to an existing file.
    '''

    path = os.environ.get('PYGECKO_MSCONVERT') or shutil.which('msconvert')
    if path and os.path.exists(path):
        return path
    return None


def msconvert(input_files: list|str, output_dir, format='mzML'):

    '''
    Takes in a list of input filenames and the output directory, converts the input files to the specified format and
    saves them in the output directory.

    Args:
        input_files (list|str): Input filename or list of input filenames.
        output_dir (str): Path to the output directory.
        format (str): Output format. Default: 'mzML'.
    '''

    output_format = format

    msconvert_path = find_msconvert()
    if msconvert_path is None:
        raise FileNotFoundError('msconvert executable not found. Set PYGECKO_MSCONVERT to its path or add '
                                'msconvert to PATH.')

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

    # Execute the msconvert command; check=True raises CalledProcessError on a failed conversion.
    subprocess.run(msconvert_cmd, check=True)
    logger.info('Conversion to %s format successful.', output_format)
