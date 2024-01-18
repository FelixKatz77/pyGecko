from setuptools import setup

setup(
    name='pyGecko',
    version='0.1.0',
    description='An open-access python libary for the parsing, processing and analysis of GC/MS and GC/FID raw data',
    url='https://github.com/shuds13/pyexample',
    author='Felix Katzenburg',
    author_email='felix.katzenburg@uni-muenster.de',
    license='MIT License',
    packages=['pygecko'],
    install_requires=[
        'brain_isotopic_distribution==1.5.14',
        'epam.indigo==1.13.0',
        'matplotlib==3.6.2',
        'numpy>=1.23.5',
        'ord_schema==0.3.37',
        'pandas==2.0.3',
        'pybaselines==1.0.0',
        'pyteomics==4.6.2',
        'rdkit==2023.3.1',
        'reportlab==4.0.5',
        'scipy==1.11.4',
        'statsmodels==0.14.0',
        'xarray==2023.6.0'
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows :: Windows 11',
        'Programming Language :: Python :: 3.10'
    ],
)