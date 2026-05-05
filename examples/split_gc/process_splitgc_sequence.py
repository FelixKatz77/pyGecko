
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from pygecko.gc_tools.injection import FID_Injection, MS_Injection
from pygecko.gc_tools.peak import FID_Peak, MS_Peak
from pygecko.parsers import SplitGC_Parser


# === SCRIPT SECTION 1: CONFIG (edit per-experiment) ===

# CONFIG: paths to the OpenLab result folders.
RSLT_PATH: Path = Path(
    r'C:\Users\flori\Doktorarbeit\22_Super_GC\data\FBS-FB-021-ALL.rslt'
)
RI_SIRSLT_PATH: Path = Path(
    r'C:\Users\flori\Doktorarbeit\22_Super_GC\data\FBS-FB-021-C7C40-Sample2.sirslt'
)

RI_PICK_PEAKS_FID_KWARGS: dict = {}
RI_PICK_PEAKS_MS_KWARGS: dict = {'prominence_ms': 50}

# CONFIG: solvent delays in minutes. The post-column splitter introduces a
# small dead-volume difference between FID and MS, so the solvent peak is
# not at exactly the same time on both detectors. The values below match
# the original test_supergc_parsers.py call; tighten or loosen if you see
# the solvent tail leaking into the integration window.
SOLVENT_DELAY_FID: float = 3.00

# CONFIG: internal standard.
IS_NAME: str = 'Trimethoxybenzene'
IS_SMILES: str = 'COC1=CC(OC)=CC(OC)=C1'
IS_RT_FID: float = 5.916  # minutes, in the FID trace
IS_RT_MS: float = 5.906   # minutes, in the MS trace
IS_RT_TOLERANCE: float = 0.05  # minutes; passed to set_internal_standard

# CONFIG: RI calibration anchor.
# ``RI_Calibration.__identify_alkanes`` walks the alkane series outward from
# this anchor, naming peaks as Cn, C(n+1), C(n-1), ... in retention-time
# order. The anchor must be a peak that actually shows up in the alkane
# standard chromatogram. Dodecane (C12) is a common, mid-range choice that
# also matches the legacy thiolation example.
#
# TODO-AFTER-FIRST-RUN: open the alkane standard's FID and MS chromatograms
# in OpenLab once, record the dodecane peak retention time on each detector,
# and update the two values below. Until then, leaving the legacy thiolation
# values as placeholders so the script structure can be validated.
RI_ANCHOR_C_COUNT: int = 12
RI_ANCHOR_RT_FID: float = 4.850
RI_ANCHOR_RT_MS: float = 4.837

# CONFIG: matching tolerances. RI tolerance follows the legacy thiolation
# script's default (20 RI units, ~0.2 min on a typical method).
RI_TOLERANCE: int = 20

# CONFIG: output. If set to a path, the long-format result table is written
# to CSV. Set to None to skip writing.
OUTPUT_CSV: Optional[Path] = Path('FBS-FB-021-ALL_yields.csv')


# === SCRIPT SECTION 2: ANALYTE TABLE (one row per well) ===
# Embedded inline so the script is self-contained. Columns mirror the Excel
# table you provided: Sample, Analyte, Analyte_SMILES. The internal standard
# is the same (TMB) for every well, so it is not repeated per row.
ANALYTE_TABLE: list[dict[str, str]] = [
    {'Sample': 'A1', 'Analyte': 'Dodecane',                            'Analyte_SMILES': 'CCCCCCCCCCCC'},
    {'Sample': 'A2', 'Analyte': 'Dodecane',                            'Analyte_SMILES': 'CCCCCCCCCCCC'},
    {'Sample': 'A3', 'Analyte': 'Dodecane',                            'Analyte_SMILES': 'CCCCCCCCCCCC'},
    {'Sample': 'B1', 'Analyte': 'Bromobenzene',                        'Analyte_SMILES': 'BrC1=CC=CC=C1'},
    {'Sample': 'B2', 'Analyte': 'Bromobenzene',                        'Analyte_SMILES': 'BrC1=CC=CC=C1'},
    {'Sample': 'B3', 'Analyte': 'Bromobenzene',                        'Analyte_SMILES': 'BrC1=CC=CC=C1'},
    {'Sample': 'C1', 'Analyte': '3,5-dimethylpyridine',                'Analyte_SMILES': 'CC1=CN=CC(C)=C1'},
    {'Sample': 'C2', 'Analyte': '3,5-dimethylpyridine',                'Analyte_SMILES': 'CC1=CN=CC(C)=C1'},
    {'Sample': 'C3', 'Analyte': '3,5-dimethylpyridine',                'Analyte_SMILES': 'CC1=CN=CC(C)=C1'},
    {'Sample': 'D1', 'Analyte': '1,2,3,4-Tetrahydroisoquinoline',      'Analyte_SMILES': 'C12=C(CNCC2)C=CC=C1'},
    {'Sample': 'D2', 'Analyte': '1,2,3,4-Tetrahydroisoquinoline',      'Analyte_SMILES': 'C12=C(CNCC2)C=CC=C1'},
    {'Sample': 'D3', 'Analyte': '1,2,3,4-Tetrahydroisoquinoline',      'Analyte_SMILES': 'C12=C(CNCC2)C=CC=C1'},
    {'Sample': 'E1', 'Analyte': 'Decanenitrile',                       'Analyte_SMILES': 'CCCCCCCCCC#N'},
    {'Sample': 'E2', 'Analyte': 'Decanenitrile',                       'Analyte_SMILES': 'CCCCCCCCCC#N'},
    {'Sample': 'E3', 'Analyte': 'Decanenitrile',                       'Analyte_SMILES': 'CCCCCCCCCC#N'},
    {'Sample': 'F1', 'Analyte': 'Mesitylene',                          'Analyte_SMILES': 'CC1=CC(C)=CC(C)=C1'},
    {'Sample': 'F2', 'Analyte': 'Mesitylene',                          'Analyte_SMILES': 'CC1=CC(C)=CC(C)=C1'},
    {'Sample': 'F3', 'Analyte': 'Mesitylene',                          'Analyte_SMILES': 'CC1=CC(C)=CC(C)=C1'},
    {'Sample': 'G1', 'Analyte': '9H-Carbazole',                        'Analyte_SMILES': 'C1(NC2=C3C=CC=C2)=C3C=CC=C1'},
    {'Sample': 'G2', 'Analyte': '9H-Carbazole',                        'Analyte_SMILES': 'C1(NC2=C3C=CC=C2)=C3C=CC=C1'},
    {'Sample': 'G3', 'Analyte': '9H-Carbazole',                        'Analyte_SMILES': 'C1(NC2=C3C=CC=C2)=C3C=CC=C1'},
    {'Sample': 'H1', 'Analyte': '1,3-Diphenylacetone',                 'Analyte_SMILES': 'O=C(CC1=CC=CC=C1)CC2=CC=CC=C2'},
    {'Sample': 'H2', 'Analyte': '1,3-Diphenylacetone',                 'Analyte_SMILES': 'O=C(CC1=CC=CC=C1)CC2=CC=CC=C2'},
    {'Sample': 'H3', 'Analyte': '1,3-Diphenylacetone',                 'Analyte_SMILES': 'O=C(CC1=CC=CC=C1)CC2=CC=CC=C2'},
    {'Sample': 'I1', 'Analyte': '4-Chloro-2-(trifluoromethyl)quinoline', 'Analyte_SMILES': 'FC(C1=NC2=CC=CC=C2C(Cl)=C1)(F)F'},
    {'Sample': 'I2', 'Analyte': '4-Chloro-2-(trifluoromethyl)quinoline', 'Analyte_SMILES': 'FC(C1=NC2=CC=CC=C2C(Cl)=C1)(F)F'},
    {'Sample': 'I3', 'Analyte': '4-Chloro-2-(trifluoromethyl)quinoline', 'Analyte_SMILES': 'FC(C1=NC2=CC=CC=C2C(Cl)=C1)(F)F'},
    {'Sample': 'J1', 'Analyte': '4-Bromo-2,2-diphenylbutyronitrile',   'Analyte_SMILES': 'N#CC(C1=CC=CC=C1)(C2=CC=CC=C2)CCBr'},
    {'Sample': 'J2', 'Analyte': '4-Bromo-2,2-diphenylbutyronitrile',   'Analyte_SMILES': 'N#CC(C1=CC=CC=C1)(C2=CC=CC=C2)CCBr'},
    {'Sample': 'J3', 'Analyte': '4-Bromo-2,2-diphenylbutyronitrile',   'Analyte_SMILES': 'N#CC(C1=CC=CC=C1)(C2=CC=CC=C2)CCBr'},
    {'Sample': 'K1', 'Analyte': 'N,N-dimethylbenzenesulfonamide',      'Analyte_SMILES': 'O=S(C1=CC=CC=C1)(N(C)C)=O'},
    {'Sample': 'K2', 'Analyte': 'N,N-dimethylbenzenesulfonamide',      'Analyte_SMILES': 'O=S(C1=CC=CC=C1)(N(C)C)=O'},
    {'Sample': 'K3', 'Analyte': 'N,N-dimethylbenzenesulfonamide',      'Analyte_SMILES': 'O=S(C1=CC=CC=C1)(N(C)C)=O'},
]


# === SCRIPT SECTION 3: PER-SAMPLE QUANTIFICATION HELPER ===
def quantify_one_sample(
    sample: str,
    analyte_name: str,
    analyte_smiles: str,
    fid_injection: FID_Injection,
    ms_injection: MS_Injection,
    ri_tolerance: int = RI_TOLERANCE,
) -> dict[str, float | str | None]:
    """Quantifies a single analyte in a single well against the TMB internal standard.

    Mirrors the matching logic from ``Analysis.__match_and_quantify`` (which
    is reaction-array-bound and therefore not directly reusable here) for the
    direct-analyte-per-well case. The flow is:

        1. Match the analyte in the MS injection by parent-ion m/z derived
           from the SMILES (``MS_Injection.match_mol``). The match also
           verifies the isotope pattern (Cl/Br compounds use M+2, others
           M+1).
        2. Take the matched MS peak's retention index and look up the
           corresponding FID peak via RI matching
           (``Injection.match_ri``). When multiple FID peaks fall inside
           the RI tolerance, ``match_ri`` automatically returns the one
           with the smallest RI deviation (excluding any peak already
           flagged as the internal standard).
        3. Quantify with ``FID_Injection.quantify(rt, method='polyarc')``
           which delegates to ``Quantification.quantify_polyarc``.

    Per the user's decision, this function does NOT fall back to direct
    retention-time matching when ``match_mol`` returns None: the result is
    NaN in that case and the status column flags the cause so the user can
    investigate parent-ion issues separately.

    Args:
        sample: Sample name (e.g. ``'A1'``); used only for the result dict.
        analyte_name: Human-readable analyte name; used only for the result.
        analyte_smiles: SMILES string of the analyte. Must produce a valid
            RDKit molecule.
        fid_injection: The FID_Injection for this sample. Must have peaks
            picked, internal standard set, and retention indices assigned.
        ms_injection: The MS_Injection for this sample. Must have peaks
            picked, internal standard set, and retention indices assigned.
        ri_tolerance: Tolerance for FID-to-MS retention-index matching, in
            RI units. Defaults to ``RI_TOLERANCE`` (20).

    Returns:
        A dict with keys ``Sample``, ``Analyte``, ``Analyte_SMILES``,
        ``RT_MS``, ``RT_FID``, ``RI``, ``Yield_pct``, ``Status``. Numeric
        fields are ``np.nan`` when matching fails. ``Status`` is one of:

          * ``'ok'`` - both MS and FID matched; yield computed.
          * ``'no_ms_match'`` - ``match_mol`` returned None (parent ion
            absent or filtered out).
          * ``'no_fid_match'`` - MS matched but no FID peak within
            ``ri_tolerance`` of the MS peak's RI.
          * ``'no_ms_ri'`` - MS peak matched but had no retention index
            assigned (shouldn't happen if ``assign_ris`` ran).
    """
    result: dict[str, float | str | None] = {
        'Sample': sample,
        'Analyte': analyte_name,
        'Analyte_SMILES': analyte_smiles,
        'RT_MS': np.nan,
        'RT_FID': np.nan,
        'RI': np.nan,
        'Yield_pct': np.nan,
        'Status': '',
    }

    # Step 1: MS parent-ion match.
    ms_match: Optional[MS_Peak] = ms_injection.match_mol(analyte_smiles)
    if ms_match is None:
        result['Status'] = 'no_ms_match'
        return result
    # match_mol can return a list when multiple candidates pass; the legacy
    # match_mol path used inside Analysis takes the lowest isotope error,
    # but in our call (default kwargs) it returns a single MS_Peak. Defend
    # against the list case anyway.
    if isinstance(ms_match, list):
        ms_match = ms_match[0]
    result['RT_MS'] = float(ms_match.rt)
    if ms_match.ri is None:
        result['Status'] = 'no_ms_ri'
        return result
    result['RI'] = float(ms_match.ri)

    # Step 2: FID retention-index match.
    fid_match: Optional[FID_Peak] = fid_injection.match_ri(
        ms_match.ri,
        tolerance=ri_tolerance,
        analyte=ms_match.analyte,
    )
    if fid_match is None:
        result['Status'] = 'no_fid_match'
        return result
    result['RT_FID'] = float(fid_match.rt)

    # Step 3: Polyarc quantification. quantify_polyarc reads peak.analyte.mol
    # for both peaks; the IS analyte was set by set_internal_standard, and
    # the analyte's analyte was set by match_ri(analyte=ms_match.analyte).
    yield_pct: int = fid_injection.quantify(fid_match.rt, method='polyarc')
    result['Yield_pct'] = float(yield_pct)
    result['Status'] = 'ok'
    return result


# === SCRIPT SECTION 4: MAIN ===
def main() -> pd.DataFrame:
    """Runs the full SplitGC analysis pipeline and returns the results DataFrame.

    Returns:
        A long-format DataFrame with one row per well in ``ANALYTE_TABLE``.
        Columns: ``Sample``, ``Analyte``, ``Analyte_SMILES``, ``RT_MS``,
        ``RT_FID``, ``RI``, ``Yield_pct``, ``Status``.
    """
    analyte_df = pd.DataFrame(ANALYTE_TABLE)
    sample_filter: set[str] = set(analyte_df['Sample'])

    # --- Step A: Load sequences ---------------------------------------------
    print(f'Loading SplitGC sequences from {RSLT_PATH} ...')
    fid_sequence, ms_sequence = SplitGC_Parser.load_sequence(
        str(RSLT_PATH),
        solvent_delay_fid=SOLVENT_DELAY_FID,
        sample_filter=sample_filter,
        pos=False,
    )
    print(
        f'  Loaded {len(fid_sequence.injections)} FID and '
        f'{len(ms_sequence.injections)} MS injections after sample_filter.'
    )

    # Sanity check: the new ANDI-Chrom CDF reader reconstructs the FID time
    # axis from actual_sampling_interval (seconds) divided by 60. If the run
    # length here looks like 600-1200 instead of 10-25, the seconds-to-
    # minutes conversion is wrong.
    first_fid: FID_Injection = next(iter(fid_sequence.injections.values()))
    fid_run_length_min: float = float(first_fid.chromatogram[0].max())
    print(
        f'  Sanity check: FID run length looks like ~{fid_run_length_min:.1f} '
        f'minutes (expected ~10-25 for typical methods).'
    )

    # --- Step B: Pick peaks --------------------------------------------------
    print('Picking peaks on FID and MS sequences ...')
    fid_sequence.pick_peaks()
    ms_sequence.pick_peaks(trace_prominence=100)

    # --- Step C: Set internal standard --------------------------------------
    print(f'Setting internal standard ({IS_NAME}) ...')
    fid_sequence.set_internal_standard(
        IS_RT_FID, tolerance=IS_RT_TOLERANCE, name=IS_NAME, smiles=IS_SMILES,
    )
    ms_sequence.set_internal_standard(
        IS_RT_MS, tolerance=IS_RT_TOLERANCE, name=IS_NAME, smiles=IS_SMILES,
    )

    # --- Step D: Load RI calibration ----------------------------------------
    print(f'Loading RI calibration from {RI_SIRSLT_PATH} ...')
    ri_calibration_fid, ri_calibration_ms = SplitGC_Parser.load_ri_calibration(
        str(RI_SIRSLT_PATH),
        c_count=RI_ANCHOR_C_COUNT,
        rt_fid=RI_ANCHOR_RT_FID,
        rt_ms=RI_ANCHOR_RT_MS,
        solvent_delay_fid=SOLVENT_DELAY_FID,
        fid_pick_peaks_kwargs=RI_PICK_PEAKS_FID_KWARGS,
        ms_pick_peaks_kwargs=RI_PICK_PEAKS_MS_KWARGS

    )

    # --- Step E: Assign retention indices -----------------------------------
    # The legacy thiolation example uses alignment=True for FID (which shifts
    # FID RTs to align the dodecane peak retention time between the
    # calibration injection and each measurement injection) and alignment=
    # False for MS. Both choices are kept here. If the SplitGC's shared
    # column makes alignment unnecessary for FID too, the failure mode is
    # silent: align_factor will simply be near zero and have no measurable
    # effect on assigned RIs.
    print('Assigning retention indices ...')
    ri_calibration_fid.assign_ris(fid_sequence, alignment=False)
    ri_calibration_ms.assign_ris(ms_sequence)

    # --- Step F: Per-well quantification ------------------------------------
    print('Matching MS->FID and quantifying yields per well ...')
    rows: list[dict[str, float | str | None]] = []
    for _, analyte_row in analyte_df.iterrows():
        sample: str = analyte_row['Sample']
        if sample not in fid_sequence.injections or sample not in ms_sequence.injections:
            rows.append({
                'Sample': sample,
                'Analyte': analyte_row['Analyte'],
                'Analyte_SMILES': analyte_row['Analyte_SMILES'],
                'RT_MS': np.nan, 'RT_FID': np.nan, 'RI': np.nan,
                'Yield_pct': np.nan, 'Status': 'sample_missing',
            })
            continue
        rows.append(quantify_one_sample(
            sample=sample,
            analyte_name=analyte_row['Analyte'],
            analyte_smiles=analyte_row['Analyte_SMILES'],
            fid_injection=fid_sequence.injections[sample],
            ms_injection=ms_sequence.injections[sample],
        ))

    results_df: pd.DataFrame = pd.DataFrame(rows)
    # Sort by well to make scrolling natural (A1, A2, A3, B1, ...).
    results_df = results_df.sort_values('Sample').reset_index(drop=True)

    # --- Step G: Report ------------------------------------------------------
    print('\n=== Results ===')
    with pd.option_context('display.max_rows', None,
                           'display.max_columns', None,
                           'display.width', 200):
        print(results_df.to_string(index=False))

    n_ok: int = int((results_df['Status'] == 'ok').sum())
    print(
        f'\nQuantified {n_ok} / {len(results_df)} samples successfully. '
        f'Failure breakdown:\n'
        f'{results_df["Status"].value_counts().to_string()}'
    )

    if OUTPUT_CSV is not None:
        results_df.to_csv(OUTPUT_CSV, index=False)
        print(f'\nResults written to {OUTPUT_CSV.resolve()}')

    return results_df


if __name__ == '__main__':
    main()