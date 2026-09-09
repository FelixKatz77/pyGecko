# pyGecko
<img src="docs/pyGecko_icon.png" alt="pyGecko_Logo" width="300" height="300"/>

> pyGecko an open-source Python library for the parsing, processing and analysis of GC-MS and GC-FID raw data.

With increasing amounts of analytical and metadata generated in HTE, data processing and analysis quickly become a 
workflow's limiting step if conducted manually. The automated processing of analytical data opens up time for chemists 
to focus on relevant outcomes, enables the standardized storage of reaction data, and facilitates the integration of 
analytical methods into closed-loop systems. Herein we present pyGecko, an open-source Python library for the parsing,
processing and analysis of GC-MS and GC-FID raw data. pyGecko offers a variety of analysis tools for the automated or 
semi-automated handling of GC measurements and sequences. This includes the interpretation of measurements in the context 
of the experiment, the automatic identification of internal standards and compound identifications based on retention 
times, the mass of a molecular ion or fragment and spectral comparison. Quantification relative to an internal standard 
can be performed for GC-FID measurements. Results of an analysis as well as chromatograms and spectra can be visualized 
and reported in standardized formats like the Open Reaction Database (ORD) schema. pyGecko is designed to be easily 
integrated into automated workflows and can be used as a stand-alone tool or as a python library.

Preprint: https://chemrxiv.org/engage/chemrxiv/article-details/66adfc465101a2ffa8001761  <br>
Paper: https://doi.org/10.1039/D4DD00347K

## Installation

> [!IMPORTANT]
> To read vendor files you need to install the msConvert tool from ProteoWizard. You can download it from [here](http://proteowizard.sourceforge.net/download.html).
> You need to specify the path to the msConvert.exe before the first run of pyGecko.

pyGecko requires Python 3.10 or newer and can be installed via pip:

```bash 
git clone https://github.com/FelixKatz77/pyGecko.git
cd pyGecko
pip install -e .
```

Optional extras: `pip install -e ".[ord]"` adds Open Reaction Database export
(`Reaction_Parser`), `".[mzml]"` adds mzML export, `".[test]"` the test dependencies and
`".[docs]"` the documentation build.

To install the exact, pinned set of dependency versions instead of the newest compatible ones,
use [uv](https://docs.astral.sh/uv/) with the committed lock file:

```bash
uv sync
```

Afterward the path to the msConvert.exe needs to be specified. This can be done by running the following command:

```bash
cd pygecko
python __init__.py
```
This will prompt you to specify the path to the msConvert.exe file:

```bash
Please provide the path to the msConvert executable or specify it in the config.ini:
```
After that pyGecko is ready to use.


## Documentation
The documentation for pyGecko can be found [here](https://pygecko.readthedocs.io/en/latest/).

## Running the tests

```bash
pip install -e ".[test,ord,mzml]"
pytest
```

Two integration tests load Agilent `.D` directories and therefore need a configured msConvert
executable; they fail without one. To skip them, run `pytest -m "not msconvert"`. The mzML export
tests need the `mzml` extra; skip them with `pytest -m "not mzml"`.

The normal offline suite includes small, attributed `.xy` and mzML excerpts from the pyGecko study.
It also checks ORD/PDF export using the corresponding plate metadata. Run it with the coverage gate:

```bash
pytest -m "not msconvert and not slow" --cov=pygecko --cov-report=term-missing --cov-fail-under=80
```

### Full study-data regressions

The complete plate tests use release 1.2 of the study dataset, pinned to
[Zenodo record 14316687](https://zenodo.org/records/14316687). Raw archives and extracted files are
kept in the ignored `.test-data/` directory and are never added to Git. Download all three
checksum-verified archives, then run the plate regressions:

```bash
python -m tests.support.fetch_zenodo
PYGECKO_REAL_DATA_DIR="$PWD/.test-data/zenodo/14316687" \
  pytest tests/real_data -m "realdata and slow"
```

These regressions process all 96 FID and all 96 mzML injections for thiolation,
Buchwald–Hartwig, and AD-HoC. Thiolation and AD-HoC reproduce the checked-in yields, retention
times, and analyte assignments exactly. Buchwald–Hartwig reproduces every assignment and retention
time; 26 of its 27 reported yields are exact. Well C9 is expected to be 67% with the current peak
overlap-border correction rather than the paper's 74%, and the test requires its `overlap` flag.
The same full run is scheduled weekly in CI and can be started with `workflow_dispatch`.

## Usage
For non-automated workflows pyGecko is best used with jupyter notebooks. The notebooks folder of the repository contains
examples for the usage of pyGecko for the quantitative analysis of reaction outcomes and spectral matching. The Python 
scripts used to perform the data processing for the publication can be found in the examples folder. GC-MS and GC-FID 
raw data for all experiments is available on Zenodo.

### Split-GC: single-injection FID + MS

For instruments that split one GC column post-column to both an MS and a Polyarc-FID detector, both traces
come from a single injection and share a retention-time axis. `SplitGC_Parser.load_sequence` reads such an
OpenLab `.rslt`/`.sirslt` folder and returns paired FID/MS sequences. Because the detectors share a time axis,
FID and MS peaks can be matched directly by **nearest retention time** (`matching='rt'` in
`Analysis.calc_plate_yield` / `Analysis.calc_plate_conv`), so no retention-index alkane standard is required;
the legacy two-machine retention-index workflow remains available via `matching='ri'` (the default). If the
result folder's `.acaml` metadata file is missing (e.g. an incomplete export), the FID injections are
enumerated directly from the `AIA/*_FID1A.cdf` files.

### Starting-material conversion and remaining starting material

In addition to product yields (`Analysis.calc_plate_yield`), pyGecko can quantify a **starting material**
relative to the internal standard:

- `Analysis.calc_plate_conv` reports **conversion** (`100 - remaining%`). By default it assumes the substrate
  was charged at the same loading as the internal standard (1 equiv); for a substrate charged in excess pass
  `equivalents` (e.g. `equivalents=1.5`) so its conversion is referenced to its actual starting amount instead
  of reading as a negative conversion. The result is floored at 0.
- `Analysis.calc_plate_rsm` reports the **remaining starting material** (the raw carbon-normalised area
  relative to the internal standard, in percent). It is reported as measured and never clamped, so an
  excess substrate can read above 100%.

See `examples/split_gc/` for a worked split-GC plate.

## Supported File Formats
pyGecko supports the following file formats:

| GC-MS         | GC-FID         |
|---------------|----------------|
| .mzML         | .xy            |
| .mzXML        | .CSV           |
| .D (Agilent)  | .cdf (ANDI/AIA)|
| .RAW (Thermo) ||
| .cdf (ANDI/AIA) ||

> [!NOTE]
> To achieve the best performance, we recommend using the .mzML file format for GC-MS data.

## Exporting Data

Processed injections and sequences can be written back out to open formats. MS data goes to mzML,
FID data to ANDI/AIA netCDF:

```python
from pygecko.parsers import (write_injection_to_mzml, write_sequence_to_mzml,
                             write_injection_to_cdf, write_sequence_to_cdf)

write_injection_to_mzml(ms_injection, 'FKB-FA-060-A1.mzML')
write_sequence_to_mzml(ms_sequence, 'exported/')       # one file per injection

write_injection_to_cdf(fid_injection, 'FBS-FA-033-A1.cdf')
write_sequence_to_cdf(fid_sequence, 'exported/')
```

mzML export needs the `mzml` extra (`pip install -e ".[mzml]"`), which pulls in
[psims](https://github.com/mobiusklein/psims). netCDF export needs nothing extra.

> [!IMPORTANT]
> An export is a record of the injection **as pyGecko holds it**, not a copy of the original
> vendor file. pyGecko's readers round m/z to nominal integer mass and keep no polarity,
> instrument or acquisition metadata, so the MS1/centroid/positive terms in the written mzML are
> the writer's defaults rather than values from the source. Data written by pyGecko reads back
> through pyGecko's own readers unchanged; it is not a faithful round-trip of the raw file.

FID data is written as netCDF rather than mzML deliberately. The PSI-MS controlled vocabulary has
no term for a flame ionization detector, and none of its chromatogram types describes one, so an
mzML export of FID data would be schema-valid but semantically wrong. ANDI/AIA (ASTM E1947/E1948)
is the chromatography standard for a detector trace, and pyGecko already reads it.

## How to Cite

If you use pyGecko in your research, please cite the following publication:

**Calibration-free quantification and automated data analysis for high-throughput reaction screening**
Felix Katzenburg, et al.
*Digital Discovery*, 2025, **4**, 384-394.
DOI: [10.1039/D4DD00347K](https://doi.org/10.1039/D4DD00347K)
