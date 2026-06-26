from pathlib import Path
from typing import Iterable, Optional

from pygecko.gc_tools import (
    FID_Injection,
    FID_Sequence,
    MS_Injection,
    MS_Sequence,
    RI_Calibration,
)
from pygecko.parsers.agilent_fid_parser import Agilent_FID_Parser
from pygecko.parsers.ms_base_parser import MS_Base_Parser


class SplitGC_Parser:
    """A coordinator parser for unified SuperGC ``.rslt`` / ``.sirslt`` folders.

    Returns parallel ``FID_Sequence`` and ``MS_Sequence`` objects keyed by
    the same sample names so they can be fed directly into the existing
    pyGecko analysis pipeline.
    """

    DEFAULT_FID_SOURCE: str = 'cdf'
    """Default ``file_source`` strategy for SuperGC data.

    Set to ``'aia_cdf'`` because OpenLab CDS always exports the FID
    chromatogram as an ANDI-Chromatography file in ``<rslt>/AIA/`` for
    SuperGC sequences. Override with ``'csv'`` to use the .CSV exports in
    the result root, or ``'auto'`` to let ``Agilent_FID_Parser`` decide.
    """

    @classmethod
    def load_sequence(
            cls,
            rslt_directory: str,
            solvent_delay_fid: float | int,
            sample_filter: Optional[Iterable[str]] = None,
            file_source: str = DEFAULT_FID_SOURCE,
            pos: bool = False,
    ) -> tuple[FID_Sequence, MS_Sequence]:
        """Loads a SuperGC result folder and returns paired FID/MS sequences.

        Args:
            rslt_directory: Path to the OpenLab ``.rslt`` folder. Must contain
                an ``AIA/`` subdirectory. An ``.acaml`` (sequence-level OpenLab
                metadata) is used when present; if it is missing (e.g. an
                incomplete export) the FID injections are enumerated directly
                from the ``AIA/*_FID1A.cdf`` filenames instead.
            solvent_delay_fid: Retention time of the solvent peak in the FID
                trace, in minutes.
            sample_filter: Optional iterable of allowed ``SampleName`` values.
                Use this to exclude cleaning/conditioning runs that OpenLab
                tags as ``SampleType='Sample'`` (and that would otherwise be
                ingested as real injections). Applied to both the FID and MS
                sequences.
            file_source: FID ingestion strategy passed through to
                ``Agilent_FID_Parser.load_sequence``. Defaults to
                ``'aia_cdf'``.
            pos: Indicates whether the sample names encode plate positions.
                Forwarded to both parsers. Defaults to False.

        Returns:
            A tuple ``(fid_sequence, ms_sequence)``.

        Raises:
            FileNotFoundError: If the ``.rslt`` folder or the ``AIA/``
                subdirectory are missing.
        """
        rslt_path = Path(rslt_directory)
        if not rslt_path.exists():
            raise FileNotFoundError(f'Result directory does not exist: {rslt_path}')
        aia_path = rslt_path / 'AIA'
        if not aia_path.exists():
            raise FileNotFoundError(
                f'AIA subdirectory not found in {rslt_path}. '
                f'A SuperGC .rslt folder must contain an AIA/ subdirectory.'
            )

        # FID side. Agilent_FID_Parser already supports sample_filter.
        fid_sequence = Agilent_FID_Parser.load_sequence(
            str(rslt_path),
            solvent_delay_fid,
            pos=pos,
            file_source=file_source,
            sample_filter=sample_filter,
        )

        ms_sequence = MS_Base_Parser.load_sequence(
            str(aia_path), pos=pos, sample_filter=sample_filter,
        )

        return fid_sequence, ms_sequence

    @classmethod
    def load_injection(
            cls,
            sirslt_directory: str,
            solvent_delay_fid: float | int,
            file_source: str = DEFAULT_FID_SOURCE,
            pos: bool = False,
    ) -> tuple[FID_Injection, MS_Injection]:
        """Loads a single-injection ``.sirslt`` folder and returns paired injections.

        ``.sirslt`` folders are structurally identical to ``.rslt`` folders
        but contain exactly one injection. Internally this delegates to
        ``load_sequence`` and unpacks the single FID/MS injection from each
        sequence.

        Args:
            sirslt_directory: Path to the OpenLab ``.sirslt`` folder.
            solvent_delay_fid: Retention time of the solvent peak in the FID
                trace, in minutes.
            solvent_delay_ms: Retention time of the solvent peak in the MS
                trace, in minutes.
            file_source: FID ingestion strategy. Defaults to ``'aia_cdf'``.
            pos: Indicates whether the sample name encodes a plate position.
                Defaults to False.

        Returns:
            A tuple ``(fid_injection, ms_injection)``.

        Raises:
            ValueError: If the folder does not contain exactly one FID and
                one MS injection.
        """
        fid_sequence, ms_sequence = cls.load_sequence(
            sirslt_directory,
            solvent_delay_fid,
            sample_filter=None,
            file_source=file_source,
            pos=pos,
        )
        if len(fid_sequence.injections) != 1 or len(ms_sequence.injections) != 1:
            raise ValueError(
                f'Expected exactly one FID and one MS injection in '
                f'{sirslt_directory}, got {len(fid_sequence.injections)} FID '
                f'and {len(ms_sequence.injections)} MS.'
            )
        fid_injection = next(iter(fid_sequence.injections.values()))
        ms_injection = next(iter(ms_sequence.injections.values()))
        return fid_injection, ms_injection

    @classmethod
    def load_ri_calibration(
            cls,
            sirslt_directory: str,
            c_count: int,
            rt_fid: float,
            rt_ms: float,
            solvent_delay_fid: float | int,
            file_source: str = DEFAULT_FID_SOURCE,
            pos: bool = False,
            fid_pick_peaks_kwargs: dict | None = None,
            ms_pick_peaks_kwargs: dict | None = None
    ) -> tuple[RI_Calibration, RI_Calibration]:
        """Loads an alkane standard ``.sirslt`` folder and returns paired RI calibrations.

        The SuperGC alkane standard is one injection acquired on a single
        column and split post-column to FID and MS. Because the two
        detectors observe (slightly) different retention times for the same
        alkane (splitter dead-volume), each detector needs its own RI
        calibration anchored to its own observed retention time.

        Args:
            sirslt_directory: Path to the OpenLab ``.sirslt`` folder
                containing the alkane standard injection.
            c_count: Carbon count of the alkane whose retention time is
                provided.
            rt_fid: Retention time of the anchor alkane in the FID trace, in
                minutes.
            rt_ms: Retention time of the anchor alkane in the MS trace, in
                minutes.
            solvent_delay_fid: Retention time of the solvent peak in the FID
                trace, in minutes.
            solvent_delay_ms: Retention time of the solvent peak in the MS
                trace, in minutes.
            file_source: FID ingestion strategy. Defaults to ``'aia_cdf'``.
            pos: Forwarded to ``load_injection``. Defaults to False.
            fid_pick_peaks_kwargs: Optional kwargs forwarded to the FID
                calibration injection's pick_peaks call. Use to override
                default peak-picking thresholds for the alkane standard.
                Defaults to None (use pyGecko defaults).
            ms_pick_peaks_kwargs: Same for MS side.


        Returns:
            A tuple ``(ri_calibration_fid, ri_calibration_ms)``.
        """
        fid_injection, ms_injection = cls.load_injection(
            sirslt_directory,
            solvent_delay_fid,
            file_source=file_source,
            pos=pos,
        )
        fid_kwargs = fid_pick_peaks_kwargs or {}
        ms_kwargs = ms_pick_peaks_kwargs or {}
        ri_calibration_fid = RI_Calibration(fid_injection, c_count, rt_fid, **fid_kwargs)
        ri_calibration_ms = RI_Calibration(ms_injection, c_count, rt_ms, **ms_kwargs)
        return ri_calibration_fid, ri_calibration_ms



if __name__ == '__main__':
    pass