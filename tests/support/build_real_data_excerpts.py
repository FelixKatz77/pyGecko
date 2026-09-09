'''Regenerate the small attributed chromatogram excerpts committed for fast tests.'''

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from pygecko.gc_tools import MS_Injection
from pygecko.parsers import MS_Base_Parser, write_injection_to_mzml


RECORD_ID = 14316687
ARCHIVE_MD5 = '32cce96527368938cae9b1a48ff3d13f'
FID_RANGE = (6.7, 7.3)
MS_RANGE = (5.8, 6.1)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def build_excerpts(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)

    fid_source = source / 'FID/FKB-FA-060-A1.xy'
    fid = np.loadtxt(fid_source, delimiter='\t')
    fid = fid[(fid[:, 0] >= FID_RANGE[0]) & (fid[:, 0] <= FID_RANGE[1])]
    fid_output = destination / 'thiolation_A1_product.xy'
    np.savetxt(fid_output, fid, delimiter='\t', fmt='%.15g')

    ms_source = source / 'MS/FKB-FA-060-A1.mzML'
    injection = MS_Base_Parser.load_injection(ms_source, pos=True)
    start, end = (value * 60_000 for value in MS_RANGE)
    scans = injection.scans.loc[start:end]
    chromatogram = np.array([scans.index / 60_000, scans.sum(axis=1)])
    excerpt = MS_Injection(
        {'SampleName': injection.sample_name}, chromatogram, peaks=None, scans=scans, pos=True)
    ms_output = destination / 'thiolation_A1_product.mzML'
    write_injection_to_mzml(excerpt, ms_output)

    metadata_outputs = {}
    for source_name, output_name in (
            ('meta_data.json', 'thiolation_meta_data.json'),
            ('plate_layout.csv', 'thiolation_plate_layout.csv')):
        output = destination / output_name
        shutil.copyfile(source / source_name, output)
        metadata_outputs[output_name] = {
            'source': f'Thiolation_Plate_GC_Data/{source_name}',
            'sha256': _sha256(output),
        }

    manifest = {
        'license': 'CC-BY-4.0',
        'record': f'https://zenodo.org/records/{RECORD_ID}',
        'archive': 'Thiolation_Plate_GC_Data.zip',
        'archive_md5': ARCHIVE_MD5,
        'files': {
            fid_output.name: {
                'source': 'Thiolation_Plate_GC_Data/FID/FKB-FA-060-A1.xy',
                'retention_time_minutes': list(FID_RANGE),
                'sha256': _sha256(fid_output),
            },
            ms_output.name: {
                'source': 'Thiolation_Plate_GC_Data/MS/FKB-FA-060-A1.mzML',
                'retention_time_minutes': list(MS_RANGE),
                'sha256': _sha256(ms_output),
            },
            **metadata_outputs,
        },
    }
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='extracted Thiolation_Plate_GC_Data directory')
    parser.add_argument(
        '--destination', type=Path, default=Path('tests/real_data/fixtures'))
    args = parser.parse_args(argv)
    build_excerpts(args.source, args.destination)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
