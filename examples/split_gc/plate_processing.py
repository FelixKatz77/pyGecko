"""Split-GC plate processing with nearest-retention-time (RT) matching.

This mirrors the legacy two-machine examples (e.g. ``thiolation_of_heteroarenes/plate_processing.py``)
but for the new **split GC**: a single injection is split post-column to BOTH an MS detector and a
Polyarc-FID detector, so the FID and MS traces share one retention-time axis. FID and MS peaks are
therefore matched by **nearest retention time** (``matching='rt'``) instead of by retention index (RI).
No alkane standard / RI calibration is required.
It uses the SAME raw data as ``process_splitgc_sequence.py`` (the ``FBS-FB-021-ALL.rslt`` SuperGC
sequence) but routes the analysis through ``Analysis.calc_plate_yield(matching='rt', ...)`` instead of a
hand-rolled per-well loop, so the matching logic lives in the core library.

Contrast with the legacy/RI example (``process_splitgc_sequence.py``): there is no ``RI_SIRSLT_PATH``
(alkane standard run), no ``RI_ANCHOR_*`` config, no ``load_ri_calibration`` step and no ``assign_ris``
step. The whole RI block disappears; only an RT-matching function and tolerance remain.
"""

import os
from pathlib import Path

from pygecko.parsers import SplitGC_Parser
from pygecko.reaction import Product_Array
from pygecko.analysis.analysis import Analysis
from pygecko.visualization.visuals import Visualization


# === CONFIG (edit per-experiment) ===

# Path to the OpenLab SuperGC result folder holding BOTH detector traces (same data as the RI example).
# The raw data is not part of the repository, so point the PYGECKO_SPLITGC_RSLT environment variable at
# your own copy, e.g.:
#   export PYGECKO_SPLITGC_RSLT=/path/to/FBS-FB-021-ALL.rslt
RSLT_PATH_ENV: str = 'PYGECKO_SPLITGC_RSLT'

# Per-well expected product SMILES as an 11x3 grid keyed A1..K3 (one analyte per row, repeated per column).
LAYOUT_CSV: Path = Path(__file__).with_name('splitgc_layout.csv')

# Solvent delay on the FID trace (minutes); trims the solvent front before peak picking.
SOLVENT_DELAY_FID: float = 3.00

# Internal standard (same for every well): Trimethoxybenzene (TMB).
IS_NAME: str = 'Trimethoxybenzene'
IS_SMILES: str = 'COC1=CC(OC)=CC(OC)=C1'
IS_RT_FID: float = 5.916  # minutes, in the FID trace
IS_RT_MS: float = 5.906   # minutes, in the MS trace
IS_RT_TOLERANCE: float = 0.05

# --- Nearest-RT matching parameters -------------------------------------------------------------------
# RT_FUNC maps an MS retention time to the EXPECTED FID retention time. Because both detectors see one
# shared column, the FID and MS retention times differ only by the small splitter dead-volume offset.
# Measured on this FBS-FB-021-ALL dataset across 27 wells the FID elutes a mean of +0.0094 min (median
# +0.010 min, range +0.002..+0.015) after the MS. A constant Analysis.constant_offset(0.010) already
# centres the match window and recovers two wells (B3, F1) whose FID peak sat just outside the +-1 s
# window at offset 0, but the offset drifts across the chromatogram, so the linear model is used instead
# (t_fid = a * t_ms + b).
RT_FUNC = Analysis.linear_drift(0.9979, +0.0232)
RT_TOLERANCE: float = 1 / 60  # half-window of the RT match, in minutes (one second)

# MS analyte-detection thresholds. pyGecko defaults are deliberately strict (max_isotopic_diff=0.055,
# min_mz_fraction=2/3). This split-GC sequence was tuned against looser values, so they are opted into
# explicitly here rather than being relaxed library-wide - they raise the rate of false-positive
# assignments and should be re-checked against known standards for any new sequence.
MAX_ISOTOPIC_DIFF: float = 0.15
MIN_MZ_FRACTION: float = 1 / 3

OUTPUT_CSV: Path = Path(__file__).with_name('FBS-FB-021-ALL_yields_rt.csv')
OUTPUT_PLATE_PNG: Path = Path(__file__).with_name('FBS-FB-021-ALL_plate_rt.png')


def main():
    '''
    Runs the split-GC plate analysis with nearest-retention-time matching and writes the results.

    Loads the paired FID/MS sequences from one SuperGC result folder, picks peaks, sets the internal
    standard on both traces, matches each well's analyte across the detectors by nearest retention time
    (matching='rt'), quantifies the yields against the internal standard (Polyarc), and writes a long-format
    CSV and a plate heatmap.

    Returns:
        np.ndarray: Structured plate array with fields 'quantity', 'rt_ms', 'rt_fid' and 'flags'.

    Raises:
        SystemExit: If the PYGECKO_SPLITGC_RSLT environment variable is not set.
    '''

    rslt_path = os.environ.get(RSLT_PATH_ENV)
    if not rslt_path:
        raise SystemExit(
            f'Set {RSLT_PATH_ENV} to the OpenLab SuperGC .rslt folder holding both detector traces, '
            f'e.g. export {RSLT_PATH_ENV}=/path/to/FBS-FB-021-ALL.rslt'
        )

    # Layout: per-well expected product SMILES (no reaction transformation / metadata needed).
    layout = Product_Array(LAYOUT_CSV)

    # Load the paired FID and MS sequences from ONE SuperGC folder. pos=True so each injection's plate
    # position is derived from its sample name (e.g. 'A1'); Analysis uses that to pair MS<->FID per well.
    sample_filter = {f'{row}{col}' for row in 'ABCDEFGHIJK' for col in '123'}
    fid_sequence, ms_sequence = SplitGC_Parser.load_sequence(
        rslt_path,
        solvent_delay_fid=SOLVENT_DELAY_FID,
        sample_filter=sample_filter,
        pos=True,
    )

    # Pick peaks (detector-appropriate kwargs per sequence).
    fid_sequence.pick_peaks()
    ms_sequence.pick_peaks(trace_prominence=100)

    # Set the internal standard on both traces (each at its own observed RT).
    fid_sequence.set_internal_standard(IS_RT_FID, tolerance=IS_RT_TOLERANCE, name=IS_NAME, smiles=IS_SMILES)
    ms_sequence.set_internal_standard(IS_RT_MS, tolerance=IS_RT_TOLERANCE, name=IS_NAME, smiles=IS_SMILES)

    # Quantify the whole plate. Per well: MS identifies the analyte peak (parent-ion m/z), the FID peak is
    # found by nearest retention time (RT_FUNC + RT_TOLERANCE), and the yield is computed against the IS
    # via the calibration-free Polyarc method. The long-format CSV is written to OUTPUT_CSV.
    yield_array = Analysis.calc_plate_yield(
        ms_sequence, fid_sequence, layout,
        path=str(OUTPUT_CSV),
        matching='rt',
        rt_func=RT_FUNC,
        rt_tolerance=RT_TOLERANCE,
        max_isotopic_diff=MAX_ISOTOPIC_DIFF,
        min_mz_fraction=MIN_MZ_FRACTION,
    )

    # Plate heatmap (labels match the 11x3 / A-K x 1-3 layout).
    Visualization.visualize_plate(
        yield_array,
        well_labels=True,
        row_labels=list('ABCDEFGHIJK'),
        col_labels=['1', '2', '3'],
        path=str(OUTPUT_PLATE_PNG),
    )

    print(f'Done. Long-format yields written to {OUTPUT_CSV}.')
    print(f'Plate heatmap written to {OUTPUT_PLATE_PNG}.')
    # Note: ORD export (Reaction_Parser.build_dataset) needs a Reaction_Array carrying a metadata JSON;
    # this analyte-only Product_Array layout has no reaction metadata, so ORD export is skipped here.

    return yield_array


if __name__ == '__main__':
    main()
