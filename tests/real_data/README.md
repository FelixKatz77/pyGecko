# Study-data regressions

The golden CSV files in `expected/` are the publication results supplied for this test work. The
small files in `fixtures/` are derived from
[pyGecko study dataset release 1.2](https://zenodo.org/records/14316687), licensed CC-BY-4.0.
`fixtures/manifest.json` records their precise source paths and checksums.

The default test suite is offline and uses only these committed files. To run all three 96-well
plates, fetch the checksum-pinned archives into the ignored local cache:

```bash
python -m tests.support.fetch_zenodo
PYGECKO_REAL_DATA_DIR="$PWD/.test-data/zenodo/14316687" \
  pytest tests/real_data -m "realdata and slow"
```

The archives and extracted raw data must remain outside Git. To regenerate the committed
thiolation excerpts after independently obtaining and extracting the source archive:

```bash
python -m tests.support.build_real_data_excerpts \
  /path/to/Thiolation_Plate_GC_Data
```
