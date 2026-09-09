'''Download the pinned pyGecko study data into an untracked local cache.'''

from __future__ import annotations

import argparse
import hashlib
import os
import stat
import tempfile
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


RECORD_ID = 14316687


@dataclass(frozen=True)
class Archive:
    filename: str
    md5: str

    @property
    def url(self) -> str:
        return (f'https://zenodo.org/api/records/{RECORD_ID}/files/'
                f'{self.filename}/content')


ARCHIVES = (
    Archive('Thiolation_Plate_GC_Data.zip', '32cce96527368938cae9b1a48ff3d13f'),
    Archive('Buchwald-Hartwig_Plate_GC_Data.zip', '12a94a07c515a91581b29fd327eadca8'),
    Archive('AD-HoC_Plate_GC_Data.zip', '75f00228ca58b4f4150ef8423d6094e1'),
)


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _retrieve(url: str, filename: Path) -> None:
    urllib.request.urlretrieve(url, str(filename))


def download_archive(
        archive: Archive,
        directory: Path,
        retrieve: Callable[[str, Path], None] = _retrieve,
) -> Path:
    '''Download one archive atomically and verify its Zenodo checksum.'''

    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / archive.filename
    if destination.exists() and _md5(destination) == archive.md5:
        return destination

    handle, temporary_name = tempfile.mkstemp(
        prefix=f'.{archive.filename}.', suffix='.part', dir=directory)
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        retrieve(archive.url, temporary)
        if _md5(temporary) != archive.md5:
            raise ValueError(f'Zenodo checksum mismatch for {archive.filename}')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def extract_archive(archive: Path, directory: Path) -> None:
    '''Extract a ZIP after rejecting traversal paths and symbolic links.'''

    directory.mkdir(parents=True, exist_ok=True)
    root = directory.resolve()
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            target = (directory / member.filename).resolve()
            try:
                target.relative_to(root)
            except ValueError as error:
                raise ValueError(f'Archive contains unsafe path: {member.filename}') from error
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise ValueError(f'Archive contains unsafe path: {member.filename}')
        bundle.extractall(directory)


def fetch_all(output: Path) -> Path:
    '''Download and extract every plate archive from the pinned record.'''

    archive_directory = output / 'archives'
    for archive in ARCHIVES:
        path = download_archive(archive, archive_directory)
        extract_archive(path, output)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output', type=Path, default=Path('.test-data/zenodo/14316687'),
        help='untracked directory in which archives and extracted data are cached')
    args = parser.parse_args(argv)
    fetch_all(args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
