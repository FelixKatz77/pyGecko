Introduction
============
pyGecko is an open-source Python library for the parsing, processing and analysis of GC-MS and GC-FID raw data. pyGecko
offers a variety of analysis tools for the automated or semi-automated handling of GC measurements and sequences. This
includes the interpretation of measurements in the context of the experiment, the automatic identification of internal
standards and compound identifications based on retention times, the mass of a molecular ion or fragment and spectral
comparison. Quantification relative to an internal standard can be performed for GC-FID measurements. Results of an
analysis as well as chromatograms and spectra can be visualized and reported in standardized formats like the Open
Reaction Database (ORD) schema. pyGecko is designed to be easily integrated into automated workflows and can be used as
a stand-alone tool or as a python library.

Installation
============
To read vendor files you need to install the msConvert tool from ProteoWizard. You can download it from [here](http://proteowizard.sourceforge.net/download.html).
You need to specify the path to the msConvert.exe before the first run of pyGecko.

pyGecko can be installed via pip:

.. code-block:: bash

   git clone https://github.com/FelixKatz77/pyGecko.git
   cd pyGecko
   pip install ./

Afterward the path to the msConvert.exe needs to be specified. This can be done by running the following command:

.. code-block:: bash

   cd pygecko
   python __init__.py

This will prompt you to specify the path to the msConvert.exe file:

.. code-block:: bash

   Please provide the path to the msConvert executable or specify it in the config.ini:

After that pyGecko is ready to use.



Usage
=====
For non-automated workflows pyGecko is best used with jupyter notebooks. The notebooks folder of the repository contains
examples for the usage of pyGecko for the quantitative analysis of reaction outcomes and spectral matching. The Python
scripts used to perform the data processing for the publication can be found in the examples folder. GC-MS and GC-FID
raw data for all experiments is available on Zenodo.

Supported file formats
======================
pyGecko supports the following file formats:

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - GC-MS
     - GC-FID
   * - .mzML
     - .xy
   * - .mzXML
     - .CSV
   * - .D (Agilent)
     -
   * - .RAW (Thermo)
     -
