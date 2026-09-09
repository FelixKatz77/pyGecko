# pyGecko Architecture

## 1. Purpose and scope

This document records the architecture of pyGecko: how the package is layered, which design choices
were made, and why. Per `CLAUDE.md` it is the place where architecture is planned and design
decisions are documented, and it is updated *after* changes land so it reflects the current design.

It has two audiences:

- **Readers**, who need to understand how raw GC data becomes a quantified reaction array.
- **Contributors**, for whom the conventions in §4, §5, §9 and the extension recipes in §7 are
  **binding**. The goal of the project is a modular library that is easy to extend and maintain, and
  that only holds if new code follows the same patterns as existing code.

This document is descriptive of the design as it stands. It is not a proposal for a redesign. Where
the current design has known problems, they are recorded in §11 rather than silently corrected here.

---

## 2. Layering

The package is a one-directional pipeline. Each layer depends only on the layers to its left.

```
  parsers/          gc_tools/          analysis/         reaction/
  ────────          ─────────          ─────────         visualization/
  vendor IO    →    domain model  →    cross-detector →  data_handling/
  + conversion      + algorithms       workflows         ─────────────
                                                         output & export
```

| Package | Responsibility |
|---|---|
| [`pygecko/parsers/`](../pygecko/parsers) | Read vendor raw data, convert it, build domain objects. The only layer that touches the filesystem or a subprocess. |
| [`pygecko/gc_tools/`](../pygecko/gc_tools) | The domain model (injections, sequences, peaks, analytes) plus all single-injection signal processing and identification. |
| [`pygecko/analysis/`](../pygecko/analysis) | Cross-detector, plate-level workflows. |
| [`pygecko/reaction/`](../pygecko/reaction) | Chemistry-side plate description: layouts, reaction SMARTS, ORD export. |
| [`pygecko/visualization/`](../pygecko/visualization), [`pygecko/data_handling/`](../pygecko/data_handling) | Rendering and reporting. |

**Rules**

- `gc_tools/` must never import from `parsers/`. Parsers construct domain objects and hand them over;
  the domain model knows nothing about file formats.
- `analysis/analysis.py` is the *only* place where an MS sequence and an FID sequence meet. Neither
  `MS_Injection` nor `FID_Injection` may reference the other detector.
- `reaction/` is chemistry-only. It does not import `gc_tools/`; `analysis/` joins the two.

There is one deliberate upward edge: domain objects expose plotting convenience methods that delegate
to `Visualization` — [`Injection.view_chromatogram`](../pygecko/gc_tools/injection/injection.py#L197),
`MS_Peak.view_mass_spectrum`. These are thin one-line delegations kept for interactive/notebook use.
Treat them as a closed set, not a licence to pull more output code into the domain model.

The edge exists at **call** time only: both methods import `Visualization` inside the method body, not
at module scope. That is what keeps it from being a true cycle — `visualization` imports `gc_tools`
downward at module scope, so a module-level import back would make `pygecko.visualization`
unimportable on its own. Keep it deferred; see §11.17.

---

## 3. Domain model

Three parallel hierarchies, each split by detector:

| Base | FID | MS |
|---|---|---|
| [`Injection`](../pygecko/gc_tools/injection/injection.py#L13) | [`FID_Injection`](../pygecko/gc_tools/injection/fid_injection.py#L13) | [`MS_Injection`](../pygecko/gc_tools/injection/ms_injection.py#L16) |
| [`GC_Sequence`](../pygecko/gc_tools/sequence/gc_sequence.py#L8) | `FID_Sequence` | `MS_Sequence` |
| [`Peak`](../pygecko/gc_tools/peak/peak.py#L6) | [`FID_Peak`](../pygecko/gc_tools/peak/fid_peak.py#L6) | [`MS_Peak`](../pygecko/gc_tools/peak/ms_peak.py#L8) |

A `GC_Sequence` holds `dict[str, Injection]` keyed by sample name; an `Injection` holds
`dict[float, Peak]` keyed by retention time. Chemical identity lives in a separate
[`Analyte`](../pygecko/gc_tools/analyte.py#L4) object attached to a peak.

### 3.1 Detector split by subclass, shared behaviour on the base

`Injection` carries everything that does not depend on the detector: peak lookup, `flag_peak`,
`match_ri`, `match_rt`, plate position, `save`. Subclasses add the detector's data representation and
its `pick_peaks` implementation.

This is what makes cross-detector work possible at all: `Analysis` identifies a compound on the MS
trace and transfers it to the FID trace through base-class methods alone. There are two
detector-independent transfer coordinates, and which one applies is a property of the *instrument*,
not of the code path:

| Coordinate | Method | Applies when | Cost |
|---|---|---|---|
| Retention index | [`match_ri`](../pygecko/gc_tools/injection/injection.py#L140) | Two separate instruments, each with its own retention-time axis | Needs an alkane ladder on both detectors |
| Retention time | [`match_rt`](../pygecko/gc_tools/injection/injection.py#L178) | One injection split post-column to both detectors (§6.H) | None — the traces already share a time axis |

> **Rule.** New behaviour that is meaningful for both detectors goes on the base class. Do not
> implement it twice in `FID_Injection` and `MS_Injection`.

### 3.2 `__slots__` on every domain class

Every domain class declares class-level type annotations *and* a matching `__slots__` tuple. A
sequence holds hundreds of injections, each carrying a full scan matrix; `__slots__` removes the
per-instance `__dict__` and prevents typo-assignment of attributes that would silently do nothing.

> **Rule.** A new attribute must be added in three places: the annotation block, `__slots__`, and
> `__init__`. Subclasses list only their *own* additional slots — see
> [`ms_injection.py:37`](../pygecko/gc_tools/injection/ms_injection.py#L37).

### 3.3 Container protocols carry the ergonomics

Sequences and injections are containers, and the API leans on that:

```python
injection = sequence['FBS-FA-033-A1']   # GC_Sequence.__getitem__ by sample name
peak      = injection[4.593]            # Injection.__getitem__ by retention time
for injection in sequence: ...          # iterates injections
len(sequence); 'sample' in sequence     # __len__, __contains__
```

See [`gc_sequence.py:35-81`](../pygecko/gc_tools/sequence/gc_sequence.py#L35-L81). Plate-oriented
access is a named method rather than an operator: `GC_Sequence.get_injection_by_pos('A1')`.

> **Rule.** Keep the dunder surface to container semantics. Anything with an argument that is not a
> key belongs in a named method.

### 3.4 Data representation and units

- **Chromatogram** — `np.ndarray` of shape `(2, N)`: row 0 is time in **minutes**, row 1 is
  intensity. Both detectors use this shape, which is why `Visualization` and `Analysis_Settings` can
  treat them uniformly (`scan_rate` is derived as `chromatogram[0,2] - chromatogram[0,1]`).
- **Scan matrix** — `MS_Injection.scans` is a `pd.DataFrame` indexed by retention time in
  **milliseconds**, with integer m/z columns and zero-filled gaps. The MS TIC is derived from it:
  `np.array([scans.index / 60000, scans.sum(axis=1)])`.
- **Heterogeneous records** are NumPy **structured arrays**, not classes:
  `MS_Peak.mass_spectrum` with fields `('mz', 'intensity', 'rel_intensity')`, `RI_Calibration.alkanes`
  with `('smiles', 'c_count', 'rt')`, and the plate result array with
  `('quantity', 'rt_ms', 'rt_fid', 'flags')`.

The mixed time units (minutes on chromatograms, milliseconds on the scans index) are a real trap.

> **Rule.** Preserve the convention rather than converting ad hoc, and state the unit in the
> docstring of every new signature that takes a time.

### 3.5 Peaks are keyed by rounded retention time

Peak dictionaries are keyed by `round(rt, 3)`. Float keys mean lookups can never rely on equality, so
all matching goes through tolerance helpers — `Utilities.check_interval(value, midpoint, tolerance)`
— and returns the closest candidate within the window (`flag_peak`, `match_ri`).

> **Rule.** Never look a peak up by computed equality. Use `flag_peak`/`match_ri`, or
> `Utilities.check_interval` if you need something new.

### 3.6 Every injection carries its processing history

`Injection.history` is an ordered list of
[`Processing_Step`](../pygecko/gc_tools/history.py) objects: the parser's load call, then every
state-changing or state-deriving operation applied since. The point is traceability — a processed
injection loaded from a `.pkl` can be traced back to the raw file and the parameters that produced
it — and a record complete enough that the same state could be reached by re-loading the raw data
and re-executing the steps. That is what makes a stored run readable as a workflow after the fact.

A step records five things: `operation` (the qualified `Class.method` name, from `func.__qualname__`),
`parameters` (the bound call arguments with defaults applied), `resolved` (see §5), `timestamp` and
the pyGecko `version`. Values are encoded to JSON-serializable form *eagerly*, at record time: a step
holding a live `Analyte` would be pickled with the injection and would silently change meaning when
that object was mutated later. A value with no JSON representation — a callable above all, such as
the `func` passed to `match_rt` — is recorded as an explicit `{'unserializable': …, 'repr': …}`
marker rather than dropped, because a missing parameter is indistinguishable from one that was never
passed. `Injection.history_to_json(path=None)` is the export.

Recording happens two ways:

- **`@records_processing`** decorates the injection methods that constitute processing:
  `set_internal_standard`, `flag_peak`, `match_ri`, `match_rt` on the base; `baseline_correction`,
  `pick_peaks`, `integrate`, `quantify` on `FID_Injection`; `pick_peaks`, `match_mz`, `match_mol` on
  `MS_Injection`. Output and accessors (`view_chromatogram`, `report`, `save`) are not processing and
  are not recorded. A call that raises records nothing.
- **`Injection.record_step(operation, parameters)`** takes provenance as plain data, for steps no
  injection method produces: the parser's load call (§6.A) and `RI_Calibration.assign_ris`, which
  mutates an injection's peaks from outside it. Passing plain data is what keeps the layering rule of
  §2 intact — `gc_tools` learns nothing about file formats or the calibration's internals.

**Nested calls record once.** `pick_peaks` calls `baseline_correction` when no processed chromatogram
exists, and `set_internal_standard` calls `flag_peak`. An `_recording` guard means only the outermost
decorated call on an injection appends a step, so the history holds no step the caller never asked
for and a replay would not execute the inner work twice. The guard suppresses *nesting*, not
*repetition*: an explicit `baseline_correction()` followed by `pick_peaks()` still records two steps.

> **Rule.** A new processing method on an injection gets `@records_processing`. A new construction
> path in a parser calls `record_step` naming a *public* parser method with arguments that method
> accepts — `tests/unit/test_injection_history.py::TestHistoryIsSufficientForReplay` binds every
> recorded parameter set against its operation's real signature, and will fail otherwise.

---

## 4. Stateless classes as algorithm namespaces

The dominant pattern in the codebase: algorithms live in classes that hold **no instance state** and
consist entirely of `@staticmethod`s (occasionally `@classmethod`s). They are namespaces, not
objects — they are never instantiated.

`Peak_Detection_FID`, `Peak_Detection_MS`, `Quantification`, `Utilities`, `Visualization`,
`Analysis`, `Reaction_Parser`, and all four parsers follow this shape.

The division of labour is consistent: **domain objects hold data and delegate; namespace classes hold
algorithms and take plain data.**

```python
# FID_Injection.pick_peaks — fid_injection.py:75
self.analysis_settings.update(**kwargs)
if not isinstance(self.processed_chromatogram, np.ndarray):
    self.baseline_correction()
peaks = Peak_Detection_FID.pick_peaks(self.processed_chromatogram, self.analysis_settings)
```

`Peak_Detection_FID.pick_peaks` receives an array and a settings object — never the `Injection`. That
keeps the algorithms independently testable and prevents cycles between the peak and injection
subpackages.

Within a namespace class, the public entry point is the only non-mangled method; its internal steps
are `__private` statics (`Peak_Detection_MS.pick_peaks` → `__detect_peaks_scipy`,
`__extract_mass_spectrum`, `__initialize_peaks`).

> **Rule.** A new algorithm is a static method on the relevant namespace class. It takes arrays plus
> an `Analysis_Settings`, not an `Injection`. Its sub-steps are name-mangled statics. Prefer adding to
> an existing namespace class over creating a new one.

---

## 5. Parameter handling: `Analysis_Settings`

[`Analysis_Settings`](../pygecko/gc_tools/analysis/analysis_settings.py#L4) is the single carrier for
every processing parameter. One instance is created per injection, in the injection's constructor,
from the chromatogram, from which it derives `scan_rate`. It carries `time_range` as a plain
setting; the scan indices for that window are derived at the call site, by
`Peak_Detection_FID.baseline_correction` against the chromatogram it is about to slice, because that
is the only place the axis being indexed is known (§11.19).

Parameters thread through the code in exactly one way:

1. **Public method takes `**kwargs`** and forwards them: `self.analysis_settings.update(**kwargs)`.
   `update` validates each key and type against the hard-coded `options` dict in `__check_settings`,
   raising `KeyError` for an unknown setting and `TypeError` for a wrong type. A typo in a keyword
   argument is therefore an error, not a silently ignored kwarg.
2. **The algorithm reads each parameter** via `settings.pop('name', computed_default)`.

Note that `Analysis_Settings.pop` does **not** remove anything ([`analysis_settings.py:88`](../pygecko/gc_tools/analysis/analysis_settings.py#L88)).
It means *"the configured value if one is set, otherwise this default"*. The name is misleading; the
behaviour is deliberate, and it is what allows defaults to be **derived from the data at the call
site** rather than fixed in the constructor:

```python
# peak_detection_ms.py:58 — default height depends on the chromatogram's own noise floor
min_height = analysis_settings.pop('height', np.min(intensities[intensities != 0]) * 50)
prominence = analysis_settings.pop('prominence_ms', 1)
```

Settings persist on the injection, so a parameter set in `baseline_correction` is still in effect in
a later `pick_peaks`.

This persistence makes it a **correctness requirement**, not a stylistic one, that a method calling
`update(**kwargs)` actually reads the settings it accepts. A method that stores a parameter and then
filters with a hard-coded value gives the caller no effect where the argument was passed *and* a
delayed effect on the next method that does read it. `MS_Injection.match_mz` had exactly this shape
until it was changed to read `min_rel_intensity` and `min_mz_fraction`.

> **Rule.** If a method calls `analysis_settings.update(**kwargs)`, every threshold it then applies
> must come from `settings.pop(...)`. Never mix stored settings and literals in one predicate.

`pop` also records what it returns, into `Analysis_Settings._resolved`. This matters because most
thresholds are never configured: the default at the `pop` call site is computed from the signal
(`prominence_fid` from the mean corrected intensity, `savgol_window` by Durbin–Watson optimisation,
`boarder_threshold` from the first difference), so the caller's kwargs are a poor record of what an
algorithm actually used. The recording decorator (§3.6) clears `_resolved` at the outermost call and
snapshots it afterwards, which is what puts real numbers in a step's `resolved` mapping. Clearing
only at the outer call is load-bearing, not incidental: it is what lets a single `pick_peaks` step
report both its own settings and those of the `baseline_correction` nested inside it.

`_resolved` is in `__slots__` (there is no `__dict__`) but is not a setting: it is absent from the
`options` dict, so `update` rejects it, and `pop` refuses any underscore-prefixed key.

> **Rule.** A new tunable parameter needs five edits in `Analysis_Settings`: the class docstring's
> attribute list, the annotation block, `__slots__`, `__init__` (initialised to `None`), and the
> `options` dict — without the last one, `update` will reject it. Put the default at the `pop` call
> site, ideally derived from the signal, not in `__init__`. `min_mz_fraction` (added for split-GC) is
> the worked example to copy.

---

## 6. The pipeline

### A. IO and conversion — `parsers/`

Format dispatch happens in exactly one place,
[`MS_Base_Parser.extract_scans_from_raw_data`](../pygecko/parsers/ms_base_parser.py#L94):

- `.mzML` → [`extract_scans_from_mzml`](../pygecko/parsers/file_readers.py) (pymzml)
- `.mzXML` → `extract_scans_from_mzxml` (pyteomics)
- `.cdf` → `extract_scans_from_cdf` (netCDF4) — ANDI/AIA open format, the second native path that
  needs no external binary. Nominal-mass binned, keeping the maximum intensity per bin.
- anything else (`.D`, `.RAW`) → [`msconvert()`](../pygecko/parsers/msconvert_wraper.py#L16) into a
  `tempfile.TemporaryDirectory`, then read back as `.mzML`

`msconvert()` is a `subprocess.run` wrapper around the external ProteoWizard executable. Its path is
read **at import time** from `pygecko/config.ini` via `configparser`
([`msconvert_wraper.py:7-10`](../pygecko/parsers/msconvert_wraper.py#L7-L10)) and is populated
interactively by running `python pygecko/__init__.py`. The path may legitimately be empty: conversion
is then unavailable, but open formats still work. This is the package's only external-binary
dependency, and it is deliberately isolated behind one function.

FID data is simpler: `FID_Base_Parser.read_xy_array` reads tab-delimited `.xy` or comma-delimited
`.CSV` with `np.loadtxt`. `Agilent_FID_Parser` additionally reads ANDI/AIA `.cdf`, reconstructing the
time axis from `actual_sampling_interval` (plus optional `actual_delay_time`) and converting to
minutes so the `(2, N)` shape of §3.4 is preserved. Which of the two it uses is chosen by the
`file_source` argument (`'csv'` for the legacy layout, `'cdf'` for split-GC exports).

Vendor parsers add **only metadata extraction** and delegate signal reading to the base parsers:
`Agilent_MS_Parser` parses `sequence.xml` inside `.D` directories, `Agilent_FID_Parser` parses Agilent
`.acaml`, both with stdlib `xml.etree.ElementTree`.

**Export is the mirror image**, in
[`file_writers.py`](../pygecko/parsers/file_writers.py): `write_injection_to_mzml` /
`write_sequence_to_mzml` and `write_injection_to_cdf` / `write_sequence_to_cdf`, module-level
functions matching the shape of `file_readers.py`. They are the first true interchange export
(§8 covers why pickle is not one). Three decisions are worth recording:

- **mzML is written with [psims](https://github.com/mobiusklein/psims), not pyteomics.** pyteomics
  is read-only for mzML — across the whole library only `fasta` and `mgf` define a `write`. psims is
  the canonical writer and validates every CV term. It is an *optional* dependency behind the `mzml`
  extra, so `_load_mzml_writer` defers the import to call time; a module-scope import would make
  `pygecko.parsers` unimportable without the extra. The deferral buys no import *time*, because
  `file_readers` imports pyteomics and pyteomics imports psims at module scope whenever it is
  installed (~0.4 s of the ~1.3 s `import pygecko.parsers`).
- **FID is written as ANDI/AIA netCDF, not mzML.** The PSI-MS CV holds no term for a flame
  ionization detector, and none of the fifteen descendants of `MS:1000626` (chromatogram type)
  describes one — they are all mass-spectrometric or instrument parameters. A chromatogram-only
  mzML *is* schema-valid (`spectrumList` is `minOccurs="0"`), so FID could be forced into one, but
  nothing in the file would identify the detector. ANDI/AIA (ASTM E1947/E1948) is the
  chromatography standard for a detector trace, and it is a format pyGecko already reads, so the
  export round-trips through `Agilent_FID_Parser.__read_cdf_file`. Because ANDI reconstructs the
  time axis from a single `actual_sampling_interval`, the writer rejects a non-uniform axis rather
  than silently distorting it.
- **The export is lossy relative to the raw file, and says so.** All three MS readers round m/z to
  nominal integer mass, and `MS_Base_Parser.initialize_injection` builds the injection from
  `{'SampleName': …}` alone — polarity, MS level, instrument, scan windows and acquisition time are
  never captured. The MS1/centroid/positive `cvParam`s are therefore writer defaults, not values
  from the source. Round-trip fidelity is claimed only against pyGecko's own readers. Intensities
  are written as float64 (the readers' own dtype) rather than the more compact float32, so a
  written file can be slightly *larger* than the vendor mzML it came from; exactness is worth more
  than the bytes. The zeros the readers insert to square off the scan matrix are dropped again on
  write, which the reader's `fillna(0)` restores.

Two loading concerns are handled at this layer rather than downstream:

- **`sample_filter`** (both `Agilent_FID_Parser.load_sequence` and `MS_Base_Parser.load_sequence`) —
  an allow-list of sample names. OpenLab tags conditioning and cleaning runs as `SampleType='Sample'`,
  so they would otherwise be ingested as real injections; filtering at load keeps that vendor quirk
  out of the domain model.
- **Missing `.acaml`** — an incomplete export has no sequence-level metadata. Rather than failing,
  `Agilent_FID_Parser` enumerates injections directly from `AIA/*_FID1A.cdf`, deriving each sample
  name from the filename and synthesising metadata with `None` fields.

### B. Signal processing — `gc_tools/peak/`

**FID** ([`peak_detection_fid.py`](../pygecko/gc_tools/peak/peak_detection_fid.py)): Savitzky–Golay
smoothing, with the window auto-tuned by a Durbin–Watson statistic (`statsmodels`), then SNIP
baseline subtraction (`pybaselines`) → `scipy.signal.find_peaks` → border detection by first-derivative
threshold → overlap resolution via `gaussian_filter1d` + `argrelmin`, which sets the `"overlap"` flag
→ Simpson integration for areas.

**MS** ([`peak_detection_ms.py`](../pygecko/gc_tools/peak/peak_detection_ms.py)): `find_peaks` on the
TIC gives candidate retention times; then every m/z trace is peak-picked independently, and a trace
peak within ±5 scans of a TIC peak contributes its intensity to that peak's mass spectrum. Relative
intensities are normalised to the base peak when the `MS_Peak` is built. MS peaks carry no baseline
correction and no area — MS is used for *identification*, FID for *quantification*.

### C. Identification — `gc_tools/analysis/`

- [`RI_Calibration`](../pygecko/gc_tools/analysis/retention_indices.py#L11) — picks peaks on an alkane
  ladder injection, seeds one known alkane by `(c_count, rt)`, walks outward assigning the rest, and
  fits `scipy.stats.linregress`. `assign_ris` accepts *either* an `Injection` or a `GC_Sequence` and
  dispatches on type. Per peak, `calculate_ri` interpolates between bracketing alkanes and falls back
  to the linear fit outside the ladder's range. `alignment=True` corrects for RT drift using the
  internal standard.
- `MS_Injection.match_mol(smiles)` — computes the molecular ion m/z with RDKit, finds candidate peaks
  containing it, and confirms with an isotope-pattern check against `brainpy.isotopic_variants`.
- `Injection.match_ri(ri, tolerance)` — the **MS → FID transfer** for two-instrument data. Retention
  index is the detector-independent coordinate that lets an MS identification be located on the FID
  trace.
- `Injection.match_rt(rt, func, tolerance)` — the **MS → FID transfer** for split-GC data (§6.H).
  `func` maps the source retention time onto the expected one in this injection's trace, absorbing
  the splitter dead-volume offset; `Analysis.constant_offset(b)` and `Analysis.linear_drift(a, b)`
  build it, and the identity mapping is used when it is omitted. Unlike `match_ri` it takes an
  `exclude_standard` flag, since the internal standard is a legitimate peak on both traces.
- [`Spectral_Match`](../pygecko/gc_tools/analysis/spectral_matching.py#L12) — library-free matching of
  two mass spectra by weighted cosine similarity (`mz**1.1 * rel_intensity**0.5`) plus an RT window.

### D. Quantification — `gc_tools/analysis/quantification.py`

All quantification is relative to the internal standard.
`quantify_polyarc` is the calibration-free default: it normalises areas by carbon count, exploiting
the FID's near-uniform per-carbon response. `quantify_calibration` uses a fitted slope/intercept.
Selection is by a `method` string in
[`FID_Injection.quantify`](../pygecko/gc_tools/injection/fid_injection.py#L111).

### E. Orchestration — `analysis/analysis.py`

[`Analysis`](../pygecko/analysis/analysis.py#L15) is the only component that sees both detectors. Per
well, `__match_and_quantify` runs:

```
plate position → expected analyte SMILES  (layout.get_product / get_substrate)
              → ms_injection.match_mol(smiles)              # identify on MS
              → matching='ri': fid_injection.match_ri(ms_peak.ri, ...)   # two instruments
                matching='rt': fid_injection.match_rt(ms_peak.rt, ...)   # split GC (§6.H)
              → __find_best_ri_match(...)                   # disambiguate candidates by
                                                            #   MS/FID height ratio to the standard
              → fid_injection.quantify(rt)                  # quantify on FID
```

The `matching` keyword selects the transfer coordinate and defaults to `'ri'`, so existing callers
are unaffected. Only the middle step differs; identification and quantification are shared.

Three quantities share this machinery, selected by an internal `mode`:

| Entry point | mode | Quantity |
|---|---|---|
| `calc_plate_yield` | `'yield'` | Product, against the internal standard |
| `calc_plate_conv` | `'conv'` | Conversion, `100 - remaining/equivalents`, floored at 0 |
| `calc_plate_rsm` | `'rsm'` | Remaining starting material, as measured and never clamped |

`equivalents` exists because the `100 % = no conversion` baseline only holds when the substrate was
charged at the same carbon-normalised amount as the standard; without it a substrate charged in
excess reads as a negative conversion. `'rsm'` is deliberately *not* clamped, so a 1.5-equiv loading
legitimately reads ~150 %.

`__find_best_ri_match` is the heart of the method: when several FID peaks fall inside the RI
tolerance, the one whose height ratio to the internal standard best matches the MS height ratio wins.

`__find_best_ri_match` is now a misnomer: it disambiguates RT candidates too. The name is kept
because renaming a name-mangled static has no functional gain (§9).

Results are returned as a **structured array**, not a bespoke result class — dtype
`[('quantity', float), ('rt_ms', float), ('rt_fid', float), ('flags', int)]`, matching the physical
well plate so it can be indexed positionally and passed straight to the heatmap. The shape is derived
from `layout.design` in `__match_and_quantify_plate`, so non-8×12 plates (the split-GC A1–K3
sequence is 11×3) work; legacy 8×12 layouts produce the same grid as before. The single-detector
`Analysis.quantify_plate` takes the same `layout` as an optional argument and falls back to 8×12
without it.

> **Rule.** Plate-level results are structured arrays with a `quantity` field and an integer `flags`
> field. Optional CSV export is a `path` keyword argument on the same method, not a separate function.

### F. Output

- [`Visualization`](../pygecko/visualization/visuals.py#L22) — `visualize_plate` (well-plate heatmap
  with optional flag markers), `view_chromatogram`, `view_mass_spectrum`, `stack_chromatograms`,
  `compare_mass_spectra` (head-to-tail). Every method takes `path=None`: it shows the figure when
  `path` is omitted and writes it when given. `visualize_plate` takes `row_labels`/`col_labels` for
  non-8×12 plates and `cbar_label` so a conversion or RSM plate is not mislabelled "Yield [%]".
- `Reaction_Parser.build_dataset` — exports to the Open Reaction Database schema (`ord_schema`
  protobufs), validated with `validations.validate_message`.
- `PDF_Report` — ReportLab document combining the heatmap, results tables, and Indigo-rendered
  structures.

#### Real-data regression layers

Study-data coverage has two layers. The default offline suite commits only narrow, attributed
chromatogram excerpts plus metadata and golden CSV results under `tests/real_data/`; the fixture
manifest records Zenodo release 1.2, the CC-BY-4.0 license, source paths, archive MD5, and excerpt
SHA-256 hashes. `tests/support/build_real_data_excerpts.py` deterministically rebuilds those files
from an extracted thiolation archive.

Whole-plate tests never access the network themselves. `tests/support/fetch_zenodo.py` separately
downloads the three pinned archives from record 14316687, verifies their published MD5 checksums,
rejects unsafe ZIP members, and extracts into the ignored `.test-data/` cache. Tests consume only
the directory named by `PYGECKO_REAL_DATA_DIR`; without it, the `realdata`/`slow` tests skip. This
keeps ordinary CI fast and offline while a weekly/manual job processes all 96 FID and 96 mzML files
per plate.

Golden comparisons are exact rather than tolerance-based because the expected CSVs were generated
by the same discrete workflow (integer yields and milliminute retention times). The sole intentional
exception is Buchwald–Hartwig C9: the current overlap-border correction produces 67% and an
`overlap` flag, whereas the publication CSV predates that correction and reports 74%.

### G. Reference usage

[`examples/buchwald_hartwig/plate_processing.py`](../examples/buchwald_hartwig/plate_processing.py) is
the canonical end-to-end script and the best single description of the intended API:

```python
rxn    = Transformation('[C,c:1][Nh1,Nh2,nh1:2].[Br,Cl:3][C,c:4]>>[C,c:1][N,n:2][C,c:4]')
layout = Reaction_Array(layout_path, rxn, meta_data_file=meta_data_path)

fid_sequence = Agilent_FID_Parser.load_sequence(fid_path, 2.7, pos=True)
ms_sequence  = MS_Base_Parser.load_sequence(ms_path, pos=True)

fid_sequence.pick_peaks()
ms_sequence.pick_peaks(prominence_ms=125)
fid_sequence.set_internal_standard(4.593, name='Dodecane', smiles='CCCCCCCCCCCC')
ms_sequence.set_internal_standard(3.324, name='Dodecane', smiles='CCCCCCCCCCCC')

ri_conf_ms.assign_ris(ms_sequence)
ri_conf_fid.assign_ris(fid_sequence, alignment=True)

yield_array = Analysis.calc_plate_yield(ms_sequence, fid_sequence, layout)
```

[`examples/split_gc/plate_processing.py`](../examples/split_gc/plate_processing.py) is the split-GC
counterpart. The shape is the same minus the RI calibration step, which the shared time axis makes
unnecessary:

```python
layout = Product_Array(LAYOUT_CSV)

fid_sequence, ms_sequence = SplitGC_Parser.load_sequence(
    rslt_path, solvent_delay_fid=3.00, sample_filter=wells, pos=True)

fid_sequence.pick_peaks()
ms_sequence.pick_peaks(trace_prominence=100)
fid_sequence.set_internal_standard(5.916, name='Trimethoxybenzene', smiles=IS_SMILES)
ms_sequence.set_internal_standard(5.906, name='Trimethoxybenzene', smiles=IS_SMILES)

yield_array = Analysis.calc_plate_yield(
    ms_sequence, fid_sequence, layout,
    matching='rt', rt_func=Analysis.linear_drift(0.9979, 0.0232), rt_tolerance=1/60)
```

> **Rule.** New workflow features must be expressible in this style: load → pick → annotate →
> analyse, with parameters passed as keyword arguments at the step where they apply.

### H. Split-GC topology

A split GC is one injection on one column, split post-column to an MS and a Polyarc-FID. Both traces
therefore originate from the same physical separation and share a retention-time axis, differing only
by a small, near-constant splitter dead-volume offset (measured at roughly +0.01 min on the reference
dataset).

This is the fact the whole `matching='rt'` path rests on. When two *separate* instruments are used,
retention times are not comparable and retention index is the only sound bridge — hence `match_ri`
and the alkane ladder. When one injection feeds both detectors, retention time is already a shared
coordinate, so the ladder is redundant and matching can be direct. The offset is absorbed by
`rt_func`; `constant_offset` suffices where it is flat, `linear_drift` where it varies across the
chromatogram.

The load path is a **coordinator parser** (§7): `SplitGC_Parser.load_sequence` reads one OpenLab
`.rslt` folder and returns *both* sequences, keyed by the same sample names so `Analysis` can pair
them per well.

> **Caution.** The two detectors derive sample identity by different routes — the MS side from the
> `.cdf` filename, the FID side from the acaml `SampleName`. They agree only when OpenLab names the
> AIA exports after the sample. `SplitGC_Parser.load_sequence` warns when the two name sets do not
> intersect, since `Analysis` can then pair no wells at all.

---

## 7. Extension points

### Adding a vendor format

Every parser satisfies the same three-method contract:

| Method | Returns |
|---|---|
| `load_sequence(directory, …)` | `MS_Sequence` / `FID_Sequence` |
| `load_injection(path, …)` | `MS_Injection` / `FID_Injection` |
| `load_ri_calibration(path, …, c_count, rt)` | `RI_Calibration` |

**This contract is duck-typed and unenforced.** There is no ABC, no `Protocol`, and no registry —
`Analysis` and the examples simply call these three names. Honour it by hand.

**Coordinator parsers are a recognised variant.**
[`SplitGC_Parser`](../pygecko/parsers/splitgc_parser.py) keeps the three method *names* but returns a
`(fid_sequence, ms_sequence)` tuple instead of a single sequence, because one `.rslt` folder holds
both detectors' data. It reads no files itself: it validates the folder layout, then delegates to
`Agilent_FID_Parser` and `MS_Base_Parser` and hands back their results. Compose existing parsers this
way rather than teaching one parser about two detectors — the single-detector parsers stay unaware of
each other, and §2's rule that only `analysis/` sees both detectors is preserved for the *domain*
objects even though a parser now loads both.

To add a vendor:

1. Write a namespace class (`<Vendor>_<Detector>_Parser`) exposing the three methods. Base parsers use
   `@staticmethod`; the Agilent parsers use `@classmethod` — either is acceptable, follow the closest
   existing parser.
2. Implement **only** the vendor's metadata extraction. Delegate signal/scan reading to
   `MS_Base_Parser.initialize_injection` or `FID_Base_Parser.read_xy_array` rather than reimplementing
   it — this is what `Agilent_FID_Parser` and `Agilent_MS_Parser` do.
3. If the format is a new *open* format rather than a vendor wrapper, add it to the dispatch in
   `MS_Base_Parser.extract_scans_from_raw_data` and to `supported_formats` in `load_sequence`.
4. Export it from [`pygecko/parsers/__init__.py`](../pygecko/parsers/__init__.py).

### Adding an export format

Exports are module-level functions in
[`file_writers.py`](../pygecko/parsers/file_writers.py), not methods on the domain objects. A
`MS_Injection.to_mzml()` convenience would need a deferred import to dodge §2's rule that
`gc_tools/` never imports `parsers/`; the functions avoid the exception entirely.

1. Add `write_injection_to_<format>(injection, path)` and, if a sequence maps onto one file per
   injection, `write_sequence_to_<format>(sequence, directory)` as a thin loop over
   `sequence.injections.values()`.
2. **Write what a pyGecko reader can read back.** The success criterion for an export is a
   round-trip through the matching `extract_scans_from_*` / `read_*` function, compared with
   `pd.testing.assert_frame_equal` or `np.testing.assert_allclose` — not a hand-checked byte
   layout. Mind §3.4's mixed units: `scans` is indexed in milliseconds, chromatograms are in
   minutes, and mzML `scan start time` is written in **minutes** because that is what
   `extract_scans_from_mzml` assumes.
3. **Keep a heavy writer library optional.** Put it behind an extra in `pyproject.toml` (see
   `mzml`), import it inside the function via a `_load_*` helper raising an `ImportError` that
   names the extra, and register a pytest marker so the tests can be deselected.
4. Do not overstate fidelity. If the domain model dropped information at read time, the docstring
   and README must say the export is a record of the injection as pyGecko holds it, not a copy of
   the source.
5. Export it from [`pygecko/parsers/__init__.py`](../pygecko/parsers/__init__.py) and add an
   `automodule` stub to `docs/source/pygecko.parsers.rst`.

### Adding a detector

Subclass `Injection`, `GC_Sequence`, and `Peak`; set `self.detector` in the injection constructor; and
implement `pick_peaks(inplace=True, **kwargs)`. `pick_peaks` is the de-facto abstract method —
`GC_Sequence.pick_peaks` calls it polymorphically on every injection, but the base `Injection` does not
declare it. Add the peak-detection algorithm as a new `Peak_Detection_<X>` namespace class.

### Adding a quantification method

Add a static method to `Quantification` taking `FID_Peak` objects, and a branch in
`FID_Injection.quantify` keyed on the `method` string.

### Adding a peak flag

`Peak.flags` is a plain `list[str]`; any flag string can be appended via `flag_peak`. A flag only needs
an entry in the [`Flags`](../pygecko/visualization/utilities.py#L5) enum if it must survive into a
plate result array — `Flags.return_flags_value` packs a flag list into the integer `flags` field that
`visualize_plate` renders. Enum values are positional, so append new members; never renumber existing
ones, or previously saved result arrays will be misread.

### Public API

Exposure is by re-export in the subpackage `__init__.py`. A class not re-exported there is internal.
Import from the subpackage (`from pygecko.parsers import Agilent_FID_Parser`), not from module paths
and not from the top-level package, which deliberately exports nothing: `pygecko/__init__.py` holds
only metadata and the interactive msConvert configuration block, so `import pygecko` stays free of the
rdkit / ord_schema / pymzml dependency tree.

---

## 8. Persistence

Injections and sequences are persisted with `pickle` / `_pickle`: `Injection.save`, `GC_Sequence.save`,
and the module-level `load_sequence` / `save_sequence`. This was chosen because a processed sequence is
a deep object graph — chromatograms, scan matrices, peaks, analytes with RDKit molecules — and pickle
round-trips it with no schema work, which suits the notebook-driven workflow the library targets.

Two consequences follow, and both are real constraints on how the code may change:

- **`.pkl` files are coupled to the class layout.** Renaming an attribute, reordering or removing a
  `__slots__` entry, or moving a class between modules breaks every previously saved sequence. This is
  the main reason the misspellings in §9 are kept.
- **Pickle executes code on load.** `.pkl` files must only be loaded from trusted sources; they are not
  an interchange format. For sharing data, use the CSV/ORD/PDF exports in §6.F. The processing history
  (§3.6) is the one structured, JSON export of an injection's own state, via `history_to_json`.

`Injection` and `Analysis_Settings` both define a `__setstate__` that fills any slot missing from the
pickled state with its default. This is what lets files written before §3.6 still load: they carry no
`history`, `_recording` or `_resolved` entry, and without the shims the first `pop` after loading such
a file raises `AttributeError`. `Analysis_Settings` needs its own because it is nested inside a pickled
injection and restores itself — `Injection.__setstate__` cannot reach it. Appending to `__slots__` is
safe for existing files; reordering or inserting is not.

---

## 9. Conventions

These are the house style. Match them in new code even where they differ from what you would write
elsewhere — internal consistency is worth more here than conformance to an external guide.

- **`Pascal_Snake_Case` class names** — `MS_Injection`, `Peak_Detection_FID`, `RI_Calibration`,
  `Analysis_Settings`. Not PEP 8, but universal in this codebase and readable for names built from
  domain acronyms. New classes match it.
- **Stateless namespace classes** for algorithms rather than bare module-level functions (§4). The
  exceptions are the pickle helpers `load_sequence` / `save_sequence` / `load_injection`, which are
  module-level by design so they can be imported without the class.
- **`__slots__` plus class-level annotations** on all domain classes (§3.2).
- **Google-style docstrings** with `Args:` / `Returns:` sections, in `'''` triple single quotes.
  `sphinx.ext.napoleon` renders them into the API docs, so every public method needs one.
- **Modern typing**: built-in generics and `X|None` unions, no `typing.Optional`. These are
  evaluated at runtime in class bodies and dataclass field annotations, so they set the package's
  hard Python 3.10 floor — see §11.16.
- **Spelling is frozen where it is public.** `boarders` (sic — borders) is the attribute name on
  `Peak` and runs through all peak-detection code; the module is `msconvert_wraper.py`. Renaming them
  would break every saved `.pkl` and every downstream script for no functional gain. Match the existing
  spelling rather than mixing both. Note the one inconsistency already present:
  [`Peak.__init__`](../pygecko/gc_tools/peak/peak.py#L32) takes the parameter as `borders` but stores it
  as `self.boarders`; subclasses pass it positionally.

---

## 10. Testing

`CLAUDE.md` sets the rules: pytest only, strict TDD (RED → GREEN → REFACTOR), ≥80% coverage overall and
100% on critical paths, external dependencies mocked, slow tests marked so `pytest -m "not slow"` stays
fast. Use the `python-testing-patterns` skill for all test work.

The suite has two halves:

- [`tests/integration/`](../tests/integration) — six modules exercising real Agilent `.D`, `.acaml`
  and `.xy` fixtures end-to-end through the parsers.
- [`tests/unit/`](../tests/unit) — one module per behaviour, plus a `conftest.py` of builders
  (`make_peak`, `make_injection`, `make_ms_injection`, `make_fid_chromatogram`,
  `make_fid_injection`, `ms_peak_factory`) that assemble domain objects from plain arrays. These
  need no fixture files and no msConvert binary. Builders live in `conftest.py` as plain module-level
  functions, imported relatively (`from .conftest import make_injection`), not as fixtures; only
  genuinely parameterised or stateful helpers are `@pytest.fixture`.

The unit half is exactly the seam the layering predicts: the namespace classes of §4 take arrays and
an `Analysis_Settings` and return data, so they can be driven with synthetic chromatograms and hand-built
mass spectra. New algorithmic code should be tested there, with an integration test only where vendor
parsing is genuinely involved.

Fixture paths in the integration half are anchored on `__file__` through
[`tests/integration/conftest.py`](../tests/integration/conftest.py), so the suite runs identically from
the repository root and from inside `tests/integration/`.

See §11.1 for the gap between these rules and the current suite.

---

## 11. Known deviations and open issues

Recorded so they are tracked rather than rediscovered. Items 1–7 are **open**: each is a statement
about the code as it stands. The subsection that follows records deviations that have since been
**resolved**, kept because the reasoning behind the fix — and, in one case, a correction to what the
defect actually did — is worth not rediscovering either.

1. **Test suite does not meet the stated coverage rules.** Coverage is now measurable and measured:
   `pytest-cov` is declared, and `pytest --cov=pygecko` reports **58%** for the suite CI runs
   against `CLAUDE.md`'s 80% requirement. The gap is concentrated in
   [`parsers/file_readers.py`](../pygecko/parsers/file_readers.py) (39%),
   [`visualization/`](../pygecko/visualization) and
   [`data_handling/reports.py`](../pygecko/data_handling/reports.py); the latter two have no direct
   tests and CI covers them only with an import smoke test. `file_readers.py` is a partial case
   since the export writers landed: the mzML reader is now exercised on every round-trip in
   `tests/integration/test_file_writers.py`, but `extract_scans_from_mzxml` and
   `extract_scans_from_cdf` remain untested because no `.mzXML` or `.cdf` fixture is checked in.
   `write_injection_to_cdf` is a way to close the `.cdf` half of that without committing a binary
   fixture.
   [`gc_tools/peak/peak_detection_fid.py`](../pygecko/gc_tools/peak/peak_detection_fid.py) left that
   list with §11.18 and is now at 94%. No `--cov-fail-under` gate is set, because
   one at 80% would keep CI permanently red — worse than no gate. Note the command in `CLAUDE.md`
   previously read `--cov=pyGecko`, which silently measured nothing (`Module pyGecko was never
   imported`); the package directory is lowercase. Two integration tests
   (`test_ms_injection.py`, `test_ms_sequence.py`) still fail wherever msConvert is absent — e.g.
   any Linux checkout — because their fixtures are Agilent `.D` directories that need conversion.
   This is deliberate: the failure is an honest signal that the binary is unconfigured, not a bug.
   They now carry a registered `msconvert` marker, so CI deselects them with
   `pytest -m "not msconvert"` while a local `pytest` behaves exactly as before. `slow` and
   `integration` markers were never actually applied to any test, so there was nothing to register;
   `--strict-markers` now makes a typo in a marker name a collection error.
2. **This document is not part of the Sphinx build.** [`docs/source/conf.py`](source/conf.py) loads only
   `autodoc`, `napoleon` and `sphinx_rtd_theme`, with no `myst_parser`, so Markdown cannot be included
   in the `toctree`. Read it directly in the repository. This is a deliberate choice: the document
   addresses contributors reading the source, not readers of the rendered API docs.
3. **`flag_peak` still keys candidates by deviation.**
   [`injection.py:110`](../pygecko/gc_tools/injection/injection.py#L110) builds
   `candidates[abs(rt - peak.rt)] = peak`, so two peaks equidistant from the target collide on one
   key. Unlike `match_ri`/`match_rt` (fixed — see §11.6 below) `flag_peak` has no `return_candidates`
   mode and always returns a single peak, so the collision only changes *which* of two equally-close
   peaks is chosen, never how many survive. Left as-is deliberately.
4. **`visuals.py` mutates global rcParams at import time.**
   [`visuals.py:187-190`](../pygecko/visualization/visuals.py#L187-L190) sets `font.family` to Arial
   process-wide when the module is imported, which affects any other plotting in the same
   interpreter and emits `findfont` warnings wherever Arial is absent — every Linux runner. The same
   module imports `pyplot` at module scope, which is why CI sets `MPLBACKEND=Agg`.
5. **`docs/source/pygecko.reaction.rst:42` autodocuments a module that no longer exists.**
   `pygecko.reaction.well_plate` was removed, but the `automodule` directive was not, so every docs
   build logs an `autodoc: failed to import` warning. Pre-existing and harmless; left for whoever
   next regenerates the `sphinx-apidoc` stubs.

6. **`time_range` is silently a no-op on `MS_Injection`.** `Peak_Detection_MS` never slices by a
   window: `pick_peaks` hands the whole chromatogram to `__detect_peaks_scipy`
   ([`peak_detection_ms.py:32`](../pygecko/gc_tools/peak/peak_detection_ms.py#L32)). `time_range` is
   still accepted, type-checked and stored by `Analysis_Settings`, so an MS caller gets no error and
   no effect. Deliberately left alone when §11.19 was fixed: `__extract_mass_spectrum`
   ([`peak_detection_ms.py:88`](../pygecko/gc_tools/peak/peak_detection_ms.py#L88)) indexes the full
   `scans` frame with the same `peak_indices` that index the chromatogram, so slicing one without
   the other desynchronises the two axes and yields a `KeyError` on `mass_spectra[rts[peak_index]]`
   or, worse, a silently mismatched mass spectrum. Windowing MS wants its own change.
7. **A cached `processed_chromatogram` makes a later `time_range` a no-op.**
   `FID_Injection.pick_peaks` only calls `baseline_correction` when `processed_chromatogram` is not
   yet an array ([`fid_injection.py:93`](../pygecko/gc_tools/injection/fid_injection.py#L93)), and
   since §11.19 that method is the only place the window is applied. So `pick_peaks()` followed by
   `pick_peaks(time_range=(5.0, 8.0))` re-picks over the *whole* chromatogram and ignores the
   window; the same holds for `savgol_window` and `max_half_window`, which also determine the
   baseline. Verified against the pre-§11.18 code, where the second call instead raised
   `IndexError: index 30347 is out of bounds for axis 0 with size 30000` — the crash was the
   double-counted offset, and removing it exposed the stale cache underneath. Not fixed with §11.19
   because the honest fix is cache invalidation — deciding which settings dirty
   `processed_chromatogram` and re-running the (expensive) baseline correction when they change —
   which is a design question, not a one-line change. Workaround: call `baseline_correction(**kwargs)`
   explicitly, or pick peaks on a freshly loaded injection.

### Resolved

6. Candidate collections in `match_ri` and `match_rt` were dictionaries keyed by absolute deviation,
   so two peaks equidistant from the target collided and one was dropped before
   `__find_best_ri_match` could weigh them by height ratio. `return_candidates=True` now returns a
   `list[Peak]` and the closest match is selected with `min(..., key=...)`;
   `__find_best_ri_match` consumes the list directly. Covered by `tests/unit/test_match_ri.py` and
   `tests/unit/test_match_rt.py`.
7. `__isotopic_ratio_check` guarded with `if not i or not j` on the index arrays returned by
   `np.where`, so a parent ion at index 0 (`array([0])`, which is falsy) made the isotope check
   report no match. It now tests `.size`. This was silently suppressing valid analyte assignments
   and was the most consequential item on this list.
8. `SplitGC_Parser.load_sequence` now warns when the FID and MS sample-name sets do not intersect.
   The two detectors still derive identity by different routes — MS from the `.cdf` filename
   (`name.split('_')[0]`), FID from the acaml `SampleName` — and that is inherent to the formats, but
   a divergence no longer surfaces as a silently empty plate. A partial overlap is legitimate (one
   detector's subset of wells) and is not reported.
9. `Analysis.quantify_plate` now takes an optional `layout` argument, deriving the grid from
   `layout.design` exactly as the MS+FID path does (§6.E). Without a layout it keeps the legacy 8×12
   grid. Note the pre-fix symptom was a `ValueError` during array assembly, not the silent
   well-dropping previously recorded here: pandas `.loc` enlargement added the out-of-grid rows and
   left the unused columns ragged.
10. `Utilities.find_empty_ranges` computed `np.isnan(signal) | np.any(signal == 0)`; `np.any`
    collapses to a scalar, so a single zero marked the entire chromatogram as empty. It now uses
    `np.isnan(signal) | (signal <= threshold)`, which also makes the documented but previously
    ignored `threshold` parameter effective. Diagnostic-only — it feeds
    `Injection._check_for_missing_signal`, not the analysis results. Because a real MS TIC carries
    30–80 genuine short dropouts per injection, that method now prints a **one-line summary** per
    injection (count, dead time, longest gap) instead of itemising every range, which would have
    put ~1700 lines on the console for a 33-well plate. Call `Utilities.find_empty_ranges` directly
    to inspect the individual gaps.
11. The dead `from xarray.util.generate_ops import inplace` in `peak_detection_ms.py` is gone, as is
    the `results='yield'` argument that `reports.py` passed to `Visualization.visualize_plate`, which
    has no such parameter.
12. Integration fixture paths were relative literals that only resolved when pytest ran from inside
    `tests/integration/`. They are now anchored on `__file__` via
    [`tests/integration/conftest.py`](../tests/integration/conftest.py), so a repo-root run — what
    `testpaths = ["tests"]` implies — behaves identically.
13. `examples/spectral_matching/spectral_matching.py` imported from the top-level package, which
    exports nothing, and raised `ImportError`; it now imports from `pygecko.parsers` like every other
    example. The top-level package deliberately still exports nothing — see §7.
14. `examples/split_gc/plate_processing.py` no longer hard-codes a personal Windows path (it reads
    `PYGECKO_SPLITGC_RSLT` and exits with a message naming the variable) and no longer assigns
    `RT_FUNC` twice. `docs/build/` and `tests/integration/.pytest_cache/` are no longer tracked.
15. **Unused and mis-declared dependencies.** `numba` was declared in
    [`pyproject.toml`](../pyproject.toml) but imported nowhere — and not even installed in the
    working venv — so it is gone. `lxml` was likewise never imported by pyGecko (both Agilent
    parsers use stdlib `xml.etree.ElementTree`), but it *is* a hard runtime import of
    `pyteomics.xml`, which `parsers/file_readers.py` imports unconditionally. Rather than drop it
    and break `.mzXML` parsing, the requirement is now expressed as `pyteomics[xml]`: pyteomics
    declares `lxml` under that extra rather than as a core dependency, so it is still installed —
    now for the true reason, and version-resolved by the package that actually needs it.
    `psycopg2-binary` was previously recorded here as **not** simply removable, because
    `ord_schema` requires `psycopg2` and dropping the binary wheel makes a fresh install compile it
    from source and fail without libpq headers. That held for `ord_schema==0.3.37`; releases from
    0.5.9 depend on `psycopg[binary,pool]>=3` instead, which ships wheels. Verified by installing
    0.5.9 in isolation: no source build, and every API surface
    [`reaction_parser.py`](../pygecko/reaction/reaction_parser.py) uses still imports with
    `UnitResolver()` instantiating. So `psycopg2-binary` is gone outright, and `ord_schema` moved to
    an `[ord]` extra floored at `>=0.5.9` — the oldest release free of both the psycopg2 build and
    the stale `protobuf<3.20` ceiling that the working venv was already violating with protobuf
    7.36.1. `Reaction_Parser` became a PEP 562 lazy module attribute in
    [`reaction/__init__.py`](../pygecko/reaction/__init__.py) so that importing anything from
    `pygecko.reaction` no longer drags in the ORD tree; `from pygecko.reaction import
    Reaction_Parser` still resolves, so the `examples/` callers were untouched. Guarded by
    [`tests/unit/test_reaction_imports.py`](../tests/unit/test_reaction_imports.py), which probes in
    a subprocess because `sys.modules` is polluted the moment any other test imports `ord_schema`.
    `netCDF4` is genuinely used, by both `.cdf` readers.
16. **Dependencies are pinned to exact versions.** All 19 runtime dependencies carried `==` pins and
    `requires-python` was `">=3.10, <3.13"`. Exact pins in a library's `install_requires` make the
    package uninstallable alongside almost anything else, and they were already fiction in practice
    — see the protobuf case in §11.15. Floors are now the previously-pinned versions, which is the
    honest claim: those are the oldest versions validated, and the ones the published work used.
    There are no upper bounds, because an upper bound in a library is a promise about software that
    does not exist yet; discovering breakage is CI's job, not the metadata's. The residual risk this
    accepts: `matplotlib`, `rdkit`, `netCDF4` and `statsmodels` declare `numpy` without an upper
    bound, so a resolve *forced* down to an old C-extension build alongside numpy ≥2 can still raise
    an ABI error. A clean resolve never picks that combination, and the two CI legs bracket it.
    `requires-python` lost its ceiling rather than gaining a `<3.14` one: it is resolver-visible, so
    a cap makes pip refuse to install pyGecko at all on a newer interpreter, whereas without one a
    failure is a legible missing-wheel error. Classifiers, not `requires-python`, are the "we test
    this" signal. The `<3.13` cap never described a code limitation — it was a consequence of
    `matplotlib==3.6.2` and `numba` having no cp312 wheels. The 3.10 floor stays, and is
    code-mandated: runtime-evaluated PEP 604 unions in class bodies and dataclass field annotations
    (§9). **No code changes were needed.** The whole modern stack was verified green before the pins
    were loosened — py3.13 with numpy 2.5.3, pandas 3.0.5, matplotlib 3.11.1, scipy 1.18.1 and rdkit
    2026.3.6 gives the same result as py3.10 on the old pins. Three breaks predicted from reading
    the upgrade notes were each disproved by running the pattern on both stacks: `set_xticklabels([])`
    after a `MultipleLocator` at [`visuals.py:262`](../pygecko/visualization/visuals.py#L262) renders
    fine; the `fillna(0, inplace=True)` at
    [`file_readers.py:35`](../pygecko/parsers/file_readers.py#L35) never sees an object column,
    because the frame is built from all-numeric dicts; and the `.loc` enlargement at
    [`analysis.py:185`](../pygecko/analysis/analysis.py#L185) indexes with a string label, not a
    tuple, and yields an identical array under pandas 3. Reproducibility moved to the layer it
    belongs on: [`uv.lock`](../uv.lock) pins a known-good environment for developers and CI, while
    `pip install -e .` resolves normally. The lock is a record going forward, not a retroactive one
    — the environment behind the publication lives in git history.
17. **`pygecko.visualization` could not be the first pyGecko import.**
    [`visuals.py:4`](../pygecko/visualization/visuals.py#L4) imports `Utilities` from
    `pygecko.gc_tools.utilities`, and importing that submodule runs the whole 20-module
    `gc_tools/__init__.py` first, which reaches `peak/ms_peak.py` — which imported `Visualization`
    back from `pygecko.visualization`, still only partway through its own first line. Hence
    `ImportError: cannot import name 'Visualization' from partially initialized module`.
    `import pygecko.gc_tools` worked only because `gc_tools.utilities` was already in `sys.modules`
    by the time `visuals.py` ran. **Six of eleven entry points failed**, two of them documented
    public API: `pygecko.data_handling` was a second, initially unrecorded break point, because
    [`reports.py:7`](../pygecko/data_handling/reports.py#L7) imports `visualization.visuals` before
    anything else, so `from pygecko.data_handling import PDF_Report` was unusable as a first import.
    [`analysis/analysis.py`](../pygecko/analysis/analysis.py) was also one edit away from breaking:
    it imports `gc_tools` on line 7 and `visualization.utilities` on line 12 and worked *only* in
    that order — swapping the two lines, exactly what an editor's "organize imports" does, broke
    `import pygecko.analysis`. Both verified, before and after.
    The fix keeps the two delegations §2 describes and defers them to call time: the
    `from pygecko.visualization import Visualization` in
    [`ms_peak.py`](../pygecko/gc_tools/peak/ms_peak.py) and
    [`injection.py`](../pygecko/gc_tools/injection/injection.py) moved inside `view_mass_spectrum`
    and `view_chromatogram`, the only places either file uses the name. Each carries a comment
    saying why, because hoisting it back to module scope silently restores the cycle. `visuals.py:4`
    deliberately stays at module scope: `visualization → gc_tools` is the correct downward direction
    per §2, and deferring *that* would have put the workaround on the layer that is supposed to
    depend downward while leaving the inverted edge in place.
    Guarded by [`tests/unit/test_import_graph.py`](../tests/unit/test_import_graph.py), which checks
    every entry point in a subprocess — in-process checks pass vacuously once `sys.modules` is warm
    — and asserts that importing `gc_tools` does not pull in `visualization`. The change is
    behaviourally neutral: both methods render byte-identical output to before.

18. **Peak boarders were reported in minutes but re-integrated as scan indices, and the time-range
    offset was counted twice.** Two defects a few lines apart in the FID detection path, both
    invisible in the default configuration.

    `Peak_Detection_FID.__detect_peaks` converted boarders from scan indices to minutes on its last
    line, *after* `__calculate_areas` had already used them as indices (which is why the areas
    `pick_peaks` sets were always right). Minutes is the intended contract:
    `Injection._check_for_peak` compares `peak.boarders` against the time values the chromatogram
    plot passes it. But [`FID_Injection.integrate`](../pygecko/gc_tools/injection/fid_injection.py)
    then did `self.chromatogram[1][round(peak.boarders[0]):round(peak.boarders[1])]` — minute values
    used as indices, so a peak at 4 min integrated `[4:4]`, an empty slice, and `simpson` raised
    `IndexError` on any real chromatogram. The method had no callers anywhere in the repository,
    examples included, which is why it had gone unnoticed. It now looks the boarders back up on the
    chromatogram's own time axis with `np.searchsorted` and integrates the **baseline-corrected**
    signal, the one `pick_peaks` integrated: all three `Quantification` methods divide one peak area
    by another, and a baseline offset does not cancel between peaks of different width, so
    integrating the raw signal would have silently shifted every yield.

    The conversion itself read
    `((peak_boarders + indices_range[0]) * scan_rate) + chrom_corr[0][0]`, and the line above it did
    `peak_indices = peak_indices + indices_range[0]`. But `chrom_corr` is the chromatogram
    `baseline_correction` already sliced by `indices_range`, so `chrom_corr[0][0]` *is* the window's
    start time and adding the offset again double-counted it; the retention times, meanwhile,
    indexed the sliced array with indices offset into the unsliced one. At the default
    `indices_range[0] == 0` both terms vanish, which is why the suite — and the golden retention
    indices in `test_fid_ri_calibration.py` — never saw it. With `time_range=(5.0, 8.0)` set, peak
    picking raised `IndexError: index 2399 is out of bounds for axis 0 with size 1199`, and any
    boarder that survived landed past the end of the run.

    Both lines now read their values straight off `chrom_corr[0]`, the slice's own time axis, rather
    than reconstructing them from `scan_rate`. That fixes the double count and makes the round trip
    `integrate` depends on **exact**. The reason the arithmetic form was not exact is *not* uneven
    sampling — the CSV fixture's 37500 points have a spacing standard deviation of `2.15e-14`, i.e.
    uniform to floating point. It is that reconstructing a value as `index * scan_rate + t0` lands
    about `3.6e-15` away from the value actually stored on the axis, which is enough for
    `np.searchsorted` to return the neighbouring scan: measured over the truncated fixture, the
    round trip was exact for only **22441 of 30000** indices, the other 7559 off by exactly +1 scan,
    shifting the narrower of its 69 peaks by up to **2.8%**. Reading off the axis round-trips
    30000/30000. It is pinned by an integration test on the real fixture (`rel=1e-12`) alongside the
    structural unit test that every boarder is a value taken from the time axis. `__find_right_boarder` may return one past the last scan — a valid slice
    bound but not a valid time — so the right boarder is clamped to the last scan.

    `indices_range` is consequently no longer read in `__detect_peaks`; its `pop` moved to
    `baseline_correction`, which is where the slice is actually taken, and where §5's recording of
    resolved parameters now picks it up (it had been read directly, bypassing `pop`).
    `Peak_Detection_MS` performs the same index-to-minute conversion arithmetically and correctly —
    it has no `indices_range` term — and is left alone: MS peaks carry no area, so nothing converts
    back. See §11.4 for the separate, still-open question of what `time_range` is relative to.

19. **`time_range` was measured from absolute zero but applied to a chromatogram truncated at the
    solvent delay.** `Analysis_Settings.__set_indices_range` converted the window with
    `convert_time_to_scan(self.time_range, scan_rate)` — indices from time zero — and
    `Peak_Detection_FID.baseline_correction` used them to slice `FID_Injection.chromatogram`, which
    the constructor had already truncated. `time_range=(5.0, 8.0)` with `solvent_delay=2.5`
    therefore analysed **7.5–10.5 min**, and could return the wrong peak rather than merely the
    wrong number of them: on the synthetic two-peak chromatogram, asking for 3.0–5.0 returned the
    peak at 5.999 instead of the one at 4.0.

    The root cause was structural. `indices_range` was *derived* state cached on
    `Analysis_Settings`, recomputed on every `update()`, and read in exactly one place — but
    `Analysis_Settings` is built from the pre-truncation array and holds no chromatogram reference,
    so the derivation could not see the axis it was indexing. Rather than thread a time origin into
    the settings object, `indices_range` was **deleted** and the conversion moved to the single
    consumer, which has the array in hand:

    ```python
    time_range = analysis_settings.pop('time_range', None)
    if time_range:
        start, end = np.searchsorted(chromatogram[0], time_range)
    else:
        start, end = 0, None
    ```

    `np.searchsorted` needs neither an origin nor `scan_rate`, and clamps a window reaching past
    either end on its own — which also closes a trap in the old code, where a `time_range` starting
    before the solvent delay produced a *negative* index that numpy silently reinterprets as
    slicing from the end. A stored origin would have needed a new slot, an `__init__` reorder (the
    old `__set_indices_range()` ran before `scan_rate` was assigned and survived only because
    `time_range` was `None` at that point), an explicit clamp, and a second `__setstate__` hook.

    Consequences worth knowing:
    - **`time_range` is now what the history records.** `resolved` carries the caller's window
      (`[5.0, 8.0]`, or `None` when unset) instead of a derived index pair — the user's intent
      rather than its translation.
    - **Removing a `__slots__` entry is exactly what §8 warns about,** so
      `Analysis_Settings.__setstate__` now skips pickled names that are no longer slots
      (`hasattr(type(self), name)`, which finds slot descriptors across the MRO). Every existing
      `.pkl` carries an `indices_range` value and would otherwise raise `AttributeError` on load;
      verified against a pickle written by the pre-change code.
    - The dead second copy of the same bug, an unused
      `Peak_Detection_FID.__set_indices_range` that reimplemented the conversion with raw float
      division and no rounding, was deleted so it cannot be wired up later.
    - `MS_Injection` needed no change, and gained nothing: see open item §11.6.
    - The stale-cache no-op this exposed is open item §11.7.

    Two latent defects in `FID_Injection.__init__` were fixed at the same time, since the change
    rewrites those lines: the truncation passed the **parameter** `solvent_delay` rather than
    `self.solvent_delay`, so the `solvent_delay=None` auto-detect branch computed a value and then
    raised `TypeError` on `None / scan_rate` — unreachable dead code, now working; and
    `if solvent_delay:` treated `solvent_delay=0` as "not provided", now `is not None`. The
    truncation arithmetic itself is unchanged, which is what keeps the 23 golden alkane retention
    times and 68 golden RIs in `test_fid_ri_calibration.py` green — they are the regression anchor
    for the default whole-chromatogram path.
