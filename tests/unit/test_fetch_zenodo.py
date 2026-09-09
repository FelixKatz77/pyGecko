import hashlib
import zipfile

import pytest

from tests.support.fetch_zenodo import (
    ARCHIVES,
    RECORD_ID,
    download_archive,
    extract_archive,
)


def test_manifest_pins_zenodo_v1_2():
    assert RECORD_ID == 14316687
    assert {archive.filename: archive.md5 for archive in ARCHIVES} == {
        'Thiolation_Plate_GC_Data.zip': '32cce96527368938cae9b1a48ff3d13f',
        'Buchwald-Hartwig_Plate_GC_Data.zip': '12a94a07c515a91581b29fd327eadca8',
        'AD-HoC_Plate_GC_Data.zip': '75f00228ca58b4f4150ef8423d6094e1',
    }


def test_download_archive_verifies_checksum_and_reuses_cache(tmp_path):
    payload = b'fixed zenodo fixture'
    archive = type(ARCHIVES[0])('fixture.zip', hashlib.md5(payload).hexdigest())
    calls = []

    def retrieve(url, filename):
        calls.append(url)
        filename.write_bytes(payload)

    path = download_archive(archive, tmp_path, retrieve=retrieve)
    cached = download_archive(archive, tmp_path, retrieve=retrieve)

    assert path == cached == tmp_path / archive.filename
    assert path.read_bytes() == payload
    assert len(calls) == 1


def test_download_archive_rejects_bad_checksum(tmp_path):
    archive = type(ARCHIVES[0])('fixture.zip', '0' * 32)

    def retrieve(url, filename):
        filename.write_bytes(b'corrupt')

    with pytest.raises(ValueError, match='checksum'):
        download_archive(archive, tmp_path, retrieve=retrieve)

    assert not (tmp_path / archive.filename).exists()


def test_extract_archive_rejects_parent_traversal(tmp_path):
    archive = tmp_path / 'unsafe.zip'
    with zipfile.ZipFile(archive, 'w') as bundle:
        bundle.writestr('../outside.txt', 'unsafe')

    with pytest.raises(ValueError, match='unsafe path'):
        extract_archive(archive, tmp_path / 'data')

    assert not (tmp_path / 'outside.txt').exists()


def test_extract_archive_extracts_safe_members(tmp_path):
    archive = tmp_path / 'safe.zip'
    with zipfile.ZipFile(archive, 'w') as bundle:
        bundle.writestr('plate/result.txt', 'ok')

    extract_archive(archive, tmp_path / 'data')

    assert (tmp_path / 'data/plate/result.txt').read_text() == 'ok'
