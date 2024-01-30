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

config = configparser.ConfigParser()
config_path = Path(__file__).parent.joinpath('config.ini')
config.read(config_path)
msconvert_exe = config.get('msConvertSettings','exe_path')
if not msconvert_exe:
    msconvert_exe = input("Please provide the path to the msConvert executable or specify it in the config.ini:")
    if pathlib.Path.exists(Path(msconvert_exe)):
        config['msConvertSettings']['exe_path'] = str(msconvert_exe)
        with open(config_path, "w") as file_object:
            config.write(file_object)
    else:
        raise FileNotFoundError("msConvert executable not found at the specified location.")
else:
    if not pathlib.Path.exists(Path(msconvert_exe)):
        raise FileNotFoundError("msConvert executable not found at the specified location. Please check the path to the msConvert executable or specify it in the config.ini.")