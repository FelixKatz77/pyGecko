"""
pyGecko

An open-source Python library for the parsing, processing and analysis of GC/MS and GC/FID raw data.
"""
import configparser
import pathlib
from pathlib import Path

__version__ = '0.1.0'
__author__ = 'Felix Katzenburg'
__credits__ = 'Muenster University'

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config_path = Path(__file__).parent.joinpath('config.ini')
    config.read(config_path)
    msconvert_exe = config.get('msConvertSettings', 'exe_path')

    # Check if the existing path is valid; if not, reset it to trigger the prompt
    if msconvert_exe and not pathlib.Path(msconvert_exe).exists():
        print(f"Current configured path not found: {msconvert_exe}")
        msconvert_exe = None

    if not msconvert_exe:
        msconvert_exe = input(
            "Please provide the path to the msConvert executable or specify it in the config.ini (Press Enter to skip): ")
        # Remove surrounding quotes if user added them
        msconvert_exe = msconvert_exe.strip('"').strip("'")

        if msconvert_exe:
            if pathlib.Path(msconvert_exe).exists():
                config['msConvertSettings']['exe_path'] = str(msconvert_exe)
                with open(config_path, "w") as file_object:
                    config.write(file_object)
            else:
                print("Warning: msConvert executable not found. Vendor file conversion will be unavailable.")
        else:
            print("Skipping msConvert configuration. Vendor file conversion will be unavailable.")