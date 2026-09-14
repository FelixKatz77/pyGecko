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

Preprint: https://chemrxiv.org/engage/chemrxiv/article-details/66adfc465101a2ffa8001761

Installation
============
To read vendor files you need to install the msConvert tool from ProteoWizard. You can download it from [here](http://proteowizard.sourceforge.net/download.html).
You need to specify the path to the msConvert.exe before the first run of pyGecko.

pyGecko is published on PyPI as ``pygecko-gc`` (the import name stays ``pygecko``):

.. code-block:: bash

   pip install pygecko-gc

Afterward tell pyGecko where msConvert is. It looks, each time a vendor file is converted, first at
the ``PYGECKO_MSCONVERT`` environment variable and then for an ``msconvert`` on your ``PATH``. Set
the variable once for your user account (on Windows via *Start → "Edit environment variables for
your account" → New*, on Linux/macOS via ``export PYGECKO_MSCONVERT=/path/to/msconvert`` in your
shell profile), or for the current session from Python:

.. code-block:: python

   import os
   os.environ["PYGECKO_MSCONVERT"] = r"C:\path\to\msconvert.exe"

After that pyGecko is ready to use. Without msConvert, open formats (``.mzML``, ``.mzXML``,
``.cdf``, ``.xy``, ``.csv``) still work.



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
