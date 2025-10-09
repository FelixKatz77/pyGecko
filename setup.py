from setuptools import setup
import setuptools

setup(
    name='pyGecko',
    version='0.1.0',
    description='An open-source Python libary for the parsing, processing and analysis of GC/MS and GC/FID raw data',
    url='https://github.com/shuds13/pyexample',
    author='Felix Katzenburg',
    author_email='felix.katzenburg@uni-muenster.de',
    license='MIT License',
    install_requires=[
        'brain_isotopic_distribution==1.5.14',
        'epam.indigo==1.13.0',
        'lxml==5.1.0',
        'matplotlib==3.6.2',
        'numba==0.58.1',
        'numpy==1.26.3',
        'ord_schema==0.3.37',
        'pandas==2.0.3',
        'pybaselines==1.0.0',
        'pymzml==2.5.2',
        'pyteomics==4.6.2',
        'rdkit==2023.3.1',
        'reportlab==4.0.5',
        'scipy==1.11.4',
        'sphinx-rtd-theme==2.0.0',
        'statsmodels==0.14.0',
        'xarray==2023.6.0',
        'scikit-learn== 1.7.2'
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
    include_package_data=True,
    packages=setuptools.find_packages(),
    package_data={'pygecko': ['config.ini']},
    python_requires='>= 3.10, < 3.11',
)