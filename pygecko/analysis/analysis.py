import numpy as np
import pandas as pd
from statistics import linear_regression

from sklearn.linear_model import LinearRegression

from pygecko.gc_tools import FID_Peak, FID_Injection
from pygecko.gc_tools.sequence import MS_Sequence, FID_Sequence
from pygecko.gc_tools.analyte import Analyte
from pygecko.reaction import Reaction_Array, Product_Array
from numpy.lib.recfunctions import unstructured_to_structured
from pygecko.visualization.utilities import Flags


class Analysis:

    '''
    A class wrapping functions to analyze GC data.
    '''

    @staticmethod
    def calc_plate_yield(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Reaction_Array|Product_Array,
                         path: str|None = None, **kwargs) -> np.ndarray:

        '''
        Matches GC-MS and GC-FID peaks and quantifies the yields of the reactions.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Reaction_Array|Product_Array): Reaction_Array or Product_Array object containing the reaction
            layout.
            path (str|None, optional): Path to write the results to. Defaults to None.
            **kwargs: Matching options forwarded to the matching routine. ``matching`` selects the MS-to-FID
                matching mode: ``'ri'`` (retention index, default; legacy two-machine workflow) or ``'rt'``
                (nearest retention time; split-GC workflow). For ``matching='rt'`` use ``rt_func`` (an MS-to-FID
                retention-time mapping, e.g. ``Analysis.constant_offset()``) and ``rt_tolerance`` (RT match
                half-window in minutes, default 1/60). For ``matching='ri'`` use ``ri_tolerance`` (RI match
                tolerance, default 20). Remaining kwargs are forwarded to ``MS_Injection.match_mol``.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''

        return Analysis.__match_and_quantify_plate(ms_sequence, fid_sequence, layout, path, mode='yield', **kwargs)

    @staticmethod
    def calc_plate_conv(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Reaction_Array,
                        path: str|None = None, index:int=0, equivalents:float=1.0, **kwargs) -> np.ndarray:

        '''
        Matches GC-MS and GC-FID peaks and quantifies the conversion of the reactions.

        Conversion is computed as ``100 - remaining%``, where ``remaining%`` is the substrate's
        carbon-normalised (Polyarc) peak area relative to the internal standard. ``remaining% = 100`` is taken
        to mean "0% conversion", which only holds if the substrate was charged at the same carbon-normalised
        molar amount as the internal standard (1 equiv). For a substrate charged at a different loading, pass
        ``equivalents`` so the baseline is scaled accordingly (``conversion = 100 - remaining% / equivalents``);
        otherwise an excess substrate reads as a negative conversion. The result is floored at 0.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Reaction_Array): Well_Plate object containing the combinatorial reaction layout.
            path (str|None, optional): Path to write the results to. Defaults to None.
            index (int, optional): Substrate index used to pick the analyte in conversion mode. Defaults to 0.
            equivalents (float, optional): Molar loading of the tracked substrate relative to the internal
                standard (which defines the 100% = no-conversion baseline). Defaults to 1.0. Pass e.g. 1.5 for
                a substrate charged at 1.5 equiv so its conversion is referenced to its actual starting amount.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''

        return Analysis.__match_and_quantify_plate(ms_sequence, fid_sequence, layout, path, mode='conv',
                                                   index=index, equivalents=equivalents, **kwargs)

    @staticmethod
    def constant_offset(b:float=0.0):

        '''
        Returns a retention-time mapping function applying a constant offset, for use as the rt_func of the
        nearest-retention-time matching mode (matching='rt').

        For a split GC the FID and MS traces originate from a single injection and share a retention-time axis
        up to a small, near-constant splitter dead-volume offset. The returned function maps a source (MS)
        retention time to the expected (FID) retention time as t_fid = t_ms + b.

        Args:
            b (float): Constant offset added to the source retention time, in minutes. Positive if the FID
                peak elutes later than the MS peak. Defaults to 0.0 (assume identical retention times).

        Returns:
            Callable[[float], float]: Function mapping a source (MS) retention time to the expected (FID)
            retention time.
        '''

        def offset(rt:float) -> float:
            return rt + b
        return offset

    @staticmethod
    def linear_drift(a:float, b:float):

        '''
        Returns a retention-time mapping function applying a linear drift, for use as the rt_func of the
        nearest-retention-time matching mode (matching='rt').

        Models a retention-time-dependent (rather than constant) difference between the detectors as
        t_fid = a * t_ms + b. Use this when the FID-vs-MS offset is observed to grow or shrink across the
        chromatogram instead of staying constant; the slope a and intercept b can be fitted from two or more
        compounds seen on both detectors (e.g. the internal standard plus a second anchor). Provided as an
        alternative to constant_offset; not used by default.

        Args:
            a (float): Slope of the linear mapping (dimensionless).
            b (float): Intercept of the linear mapping, in minutes.

        Returns:
            Callable[[float], float]: Function mapping a source (MS) retention time to the expected (FID)
            retention time.
        '''

        def drift(rt:float) -> float:
            return a * rt + b
        return drift

    @staticmethod
    def __match_and_quantify_plate(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Reaction_Array,
                                   path: str|None = None, mode:str='yield', index:int=0, **kwargs) -> np.ndarray:

        '''
        Matches GC-MS and GC-FID peaks and quantifies the yields/conversion of the reactions.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Reaction_Array): Well_Plate object containing the combinatorial reaction layout.
            path (str|None, optional): Path to write the results to. Defaults to None.
            mode (str, optional): Parameter (yield or conversion) to quantify. Defaults to 'yield'.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''



        results_dict = Analysis.__match_and_quantify(ms_sequence, fid_sequence, layout, mode, index, **kwargs)

        # Plate shape is taken from the layout (rows A.. , columns 1..) so non-8x12 plates (e.g. the split-GC
        # A1..K3 sequence) are handled. Legacy 8x12 layouts yield the same A-H / 1-12 grid as before.
        row_labels = [str(x) for x in layout.design.x.values]
        col_labels = [str(y) for y in layout.design.y.values]
        results_df = pd.DataFrame(columns=col_labels, index=row_labels)

        for key, value in results_dict.items():
            results_df.loc[key[0], key[1:]] = value

        results_array = results_df.to_numpy()
        results_array = np.array([row.tolist() for row in results_array])
        dtype = np.dtype([('quantity', float), ('rt_ms', float), ('rt_fid', float), ('flags', int)])
        results_array = unstructured_to_structured(results_array[:,:,:-1], dtype=dtype)
        if path:
            if mode == 'yield':
                quantity = 'Yield [%]'
                report_df = pd.DataFrame.from_dict(results_dict, orient='index',
                                                   columns=[quantity, 'RT-MS [min]', 'RT-FID [min]', 'Flags', 'Analyte'])
                report_df.drop(columns=['Flags'], inplace=True)
            else:
                quantity = 'Conversion [%]'
                report_df = pd.DataFrame.from_dict(results_dict, orient='index',
                                                   columns=[quantity, 'RT-MS [min]', 'RT-FID [min]', 'Flags', 'Analyte'])
                report_df.drop(columns=['Flags'], inplace=True)
            report_df.sort_index(inplace=True)
            report_df.to_csv(path)
        return results_array

    @staticmethod
    def __match_and_quantify(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Reaction_Array, mode: str,
                             index: int = 0, **kwargs) -> dict[str, list[float, str]]:

        '''
        Matches GC-MS and GC-FID peaks and quantifies the yields/conversion of the reactions for each well.

        For every MS injection the expected analyte (from the layout) is identified in the MS trace by
        parent-ion m/z (match_mol), the corresponding FID peak is located, and the yield is quantified against
        the internal standard. The MS-to-FID matching is selected by the ``matching`` keyword: ``'ri'`` bridges
        the detectors via retention index (legacy two-machine data), ``'rt'`` matches by nearest retention time
        (split-GC data, where both traces come from one injection).

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Reaction_Array): Reaction_Array or Product_Array object containing the reaction layout.
            mode (str): Parameter to quantify, 'yield' or 'conv'.
            index (int): Substrate index used to pick the analyte in conversion mode. Defaults to 0.
            **kwargs: matching (str): MS-to-FID matching mode, 'ri' (retention index, default) or 'rt' (nearest
                retention time). ri_tolerance (int): RI match tolerance for matching='ri'. Defaults to 20.
                rt_func (callable|None): MS-to-FID retention-time mapping for matching='rt' (e.g.
                Analysis.constant_offset()); the identity mapping is used if None. rt_tolerance (float): RT
                match half-window in minutes for matching='rt'. Defaults to 1/60 (one second). Remaining kwargs
                are forwarded to MS_Injection.match_mol.

        Returns:
            dict: Dictionary containing the yields/conversion, retention times, flags and the analyte smiles for
            each well, keyed by plate position.
        '''

        matching = kwargs.pop('matching', 'ri')
        ri_tolerance = kwargs.pop('ri_tolerance', 20)
        rt_func = kwargs.pop('rt_func', None)
        rt_tolerance = kwargs.pop('rt_tolerance', 1 / 60)
        equivalents = kwargs.pop('equivalents', 1.0)

        results_dict = {}
        for ms_injection in ms_sequence:
            pos = ms_injection.get_plate_position()
            fid_injection = fid_sequence.get_injection_by_pos(pos)

            # An injection without an internal standard (e.g. the IS peak was missing and a warning was raised
            # by set_internal_standard) cannot be quantified; report it as unmatched instead of crashing.
            if ms_injection.internal_standard is None or fid_injection.internal_standard is None:
                results_dict[pos] = [np.nan, np.nan, np.nan, 0, '']
                continue

            if mode == 'yield':
                analyte = layout.get_product(ms_injection.get_plate_position())
                pass
            if mode == 'conv':
                analyte = layout.get_substrate(pos, index=index)
            mz_match = ms_injection.match_mol(analyte, **kwargs)
            if mz_match:
                ms_height_ratio = mz_match.height / ms_injection[ms_injection.internal_standard.rt].height
                if matching == 'rt':
                    # Split-GC: FID and MS share one retention-time axis, so match the FID peak by nearest
                    # retention time (rt_func maps the MS rt to the expected FID rt; identity by default).
                    fid_candidates = fid_injection.match_rt(mz_match.rt, func=rt_func, tolerance=rt_tolerance,
                                                            analyte=mz_match.analyte, return_candidates=True)
                else:
                    # Legacy two-machine: bridge the detectors via retention index.
                    fid_candidates = fid_injection.match_ri(mz_match.ri, analyte=mz_match.analyte,
                                                            return_candidates=True, tolerance=ri_tolerance)

                if fid_candidates:
                    fid_match = Analysis.__find_best_ri_match(fid_candidates, fid_injection, ms_height_ratio, analyte=mz_match.analyte)
                    yield_ = fid_injection.quantify(fid_match.rt)
                    if mode == 'conv':
                        # Reference the remaining substrate to its charged loading (equivalents) so an excess
                        # substrate is not reported as negative conversion, and floor at 0.
                        yield_ = max(0.0, round(100 - yield_ / equivalents, 0))

                    flags = list(set(mz_match.flags + fid_match.flags))
                    flags = Flags.return_flags_value(flags)
                    results_dict[pos] = [yield_, mz_match.rt, fid_match.rt, flags, analyte]
                else:
                    results_dict[pos] = [np.nan, np.nan, np.nan, 0, '']
            else:
                results_dict[pos] = [np.nan, np.nan, np.nan, 0, '']
        return results_dict


    @staticmethod
    def quantify_analyte(fid_sequence:FID_Sequence, rt_analyte:float, analyte:Analyte|None=None, method='polyarc', path: str|None = None, **kwargs) -> dict:
        '''
        Returns the yields of Analytes at a given retention time in a FID sequence.

        Args:
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            rt (float): Retention time to quantify the Analytes at.
            analyte (Analyte): Analyte to quantify.
            method (str): Method to use for quantification. Defaults to 'polyarc' can be set to calibration to use calibration curve.
            path (str|None): Path to write the results to. Defaults to None.

        Returns:
            dict: Dictionary containing the yields of the Analytes in the FID sequence.
        '''



        yield_dict = {}
        for injection in fid_sequence.injections.values():
            peak = injection.flag_peak(rt_analyte, analyte=analyte)
            if peak:
                yield_ = injection.quantify(peak.rt, method=method, **kwargs)
                flags = Flags.return_flags_value(peak.flags)
                rt = peak.rt
                if peak.analyte:
                    analyte_name = peak.analyte.name
                else:
                    analyte_name = ''
            else:
                yield_ = np.nan
                flags = 0
                rt = np.nan
                analyte_name = ''
            if injection.plate_pos:
                yield_dict[injection.plate_pos] = [yield_, rt, flags, analyte_name]
            else:
                yield_dict[injection.sample_name] = [yield_, rt, flags, analyte_name]

        if path:
            report_df = pd.DataFrame.from_dict(yield_dict, orient='index',
                                               columns=['Yield [%]', 'RT-FID [min]', 'Flags', 'Analyte'])
            report_df.drop(columns=['Flags'], inplace=True)
            report_df.sort_index(inplace=True)
            report_df.to_csv(path)

        return yield_dict

    @staticmethod
    def __find_best_ri_match(match_candidates:dict[float:FID_Peak], fid_injection:FID_Injection, ms_height_ratio:float, analyte:str|None=None) -> FID_Peak:

        '''
        Returns the best retention index match based on the height ratios of the peaks.

        Args:
            ri_match (FID_Peak): FID_Peak object to match.
            fid_injection (FID_Injection): FID_Injection object to match the peak to.
            ms_height_ratio (float): Height ratio of the MS peak.
            analyte(str|None): Analyte object to be assigned to the peak. Default is None.

        Returns:
            FID_Peak: FID_Peak object with the best match based on the height ratios.
        '''

        ratio_candidates = {}
        for ri_diff, peak in match_candidates.items():
            fid_height_ratio = peak.height / fid_injection[fid_injection.internal_standard.rt].height
            height_ratio_diff = abs(fid_height_ratio  - ms_height_ratio)
            ratio_candidates[height_ratio_diff] = peak
        best_match = ratio_candidates[min(ratio_candidates)]
        if analyte:
            best_match.analyte = analyte
        return best_match

    @staticmethod
    def fit_calibration_curve(
            fid_sequence: FID_Sequence,
            analyte: Analyte,
            ratios: list[float],
            intercept: bool = False,  # set to False to force through (0,0)
    ) -> tuple[float, float]:
        """
        Fits a calibration curve to the analyte in a FID sequence using scikit-learn.

        Args:
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            analyte (Analyte): Analyte to fit the calibration curve to.
            ratios (list[float]): List of concentration ratios for calibration.
            intercept (bool): If True, fit intercept (default True). If False, force through origin.

        Returns:
            tuple[float, float]: (slope, intercept) of calibration curve.
        """
        area_ratios = []
        for injection in fid_sequence.injections.values():
            analyte_peak = injection.flag_peak(rt=analyte.rt, analyte=analyte)
            standard_peak = injection.flag_peak(injection.internal_standard.rt)
            if analyte_peak:
                area_ratio = analyte_peak.area / standard_peak.area
                area_ratios.append(area_ratio)
            else:
                area_ratios.append(0.0)

        # Prepare data for sklearn (X must be 2D)
        X = np.asarray(area_ratios, dtype=float).reshape(-1, 1)
        y = np.asarray(ratios, dtype=float)

        # Fit linear model (optionally without intercept to force through origin)
        model = LinearRegression(fit_intercept=intercept)
        model.fit(X, y)

        slope = float(model.coef_[0])
        b = float(model.intercept_) if intercept else 0.0

        # R^2 score
        r2 = float(model.score(X, y))

        print(f'Calibration fitted: Slope: {slope}, Intercept: {b}, R^2: {r2}')
        if r2 < 0.9:
            print('Calibration curve fit is not accurate.')

        return slope, b

    @staticmethod
    def quantify_plate(fid_sequence:FID_Sequence, rt_analyte:float, analyte:Analyte|None=None, method='polyarc', path: str|None = None, **kwargs) -> np.ndarray:
        '''
        Returns the yields of Analytes at a given retention time in the FID sequence of a plate.

        Args:
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            rt (float): Retention time to quantify the Analytes at.
            analyte (Analyte): Analyte to quantify.
            method (str): Method to use for quantification. Defaults to 'polyarc' can be set to calibration to use calibration curve.
            path (str|None): Path to write the results to. Defaults to None.

        Returns:
            np.ndarray: Array containing the yields of the Analytes in the FID sequence.
        '''



        yield_dict = {}
        for injection in fid_sequence.injections.values():
            peak = injection.flag_peak(rt_analyte, analyte=analyte)
            if peak:
                yield_ = injection.quantify(peak.rt, method=method, **kwargs)
                flags = Flags.return_flags_value(peak.flags)
                rt = peak.rt
                if peak.analyte:
                    analyte_name = peak.analyte.name
                else:
                    analyte_name = ''
            else:
                yield_ = np.nan
                flags = 0
                rt = np.nan
                analyte_name = ''
            if injection.plate_pos:
                yield_dict[injection.plate_pos] = [yield_, rt, flags, analyte_name]
            else:
                yield_dict[injection.sample_name] = [yield_, rt, flags, analyte_name]

        results_df = pd.DataFrame(columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
                                  index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

        for key, value in yield_dict.items():
            results_df.loc[key[0], key[1:]] = value

        results_array = results_df.to_numpy()
        results_array = np.array([row.tolist() for row in results_array])
        dtype = np.dtype([('quantity', float), ('rt_fid', float), ('flags', int)])
        results_array = unstructured_to_structured(results_array[:, :, :-1], dtype=dtype)

        if path:
            report_df = pd.DataFrame.from_dict(yield_dict, orient='index',
                                               columns=['Yield [%]', 'RT-FID [min]', 'Flags', 'Analyte'])
            report_df.drop(columns=['Flags'], inplace=True)
            report_df.sort_index(inplace=True)
            report_df.to_csv(path)

        return results_array

    @staticmethod
    def match_mz_plate(ms_sequence: MS_Sequence, layout: Reaction_Array, path: str|None = None,
                       return_canidates:bool = True, check_iso:bool = False, **kwargs) -> pd.DataFrame:

        '''
        Matches m/z values of the MS peaks to the analytes in the reaction layout.
        Args:
            ms_sequence: MS_Sequence object containing the GC-MS data.
            layout: Reaction_Array object containing the combinatorial reaction layout.
            path: (str|None): Path to write the results to. Defaults to None.
            **kwargs: Additional keyword arguments to pass to the match_mol method.

        Returns:
            pd.DataFrame: DataFrame containing the m/z values, retention times and smiles for the analytes.
        '''

        results_dict = Analysis.__match_mz(ms_sequence, layout, return_canidates=return_canidates, check_iso=check_iso,
                                           **kwargs)
        data = {}
        columns = ['RT-MS [min]']*max([len(entry[0]) for entry in results_dict.values()])
        columns.append('Analyte')
        for pos, mz_match in results_dict.items():
            if mz_match[0]:
                a = list(mz_match[0].keys())
                row = list(mz_match[0].keys()) + [np.nan]*(len(columns)-len(mz_match[0])-1)
                row.append(mz_match[1])
                data[pos] = row
            else:
                data[pos] = [np.nan]*len(columns)
        report_df = pd.DataFrame.from_dict(data, orient='index',
                                                   columns=columns)
        report_df.sort_index(inplace=True)
        if path:
            report_df.to_csv(path)
        return report_df

    @staticmethod
    def __match_mz(ms_sequence: MS_Sequence, layout: Reaction_Array, return_canidates:bool = True, check_iso:bool = False, **kwargs) -> dict[str, list[float, str]]:

        '''
        Matches m/z values of the MS peaks to the analytes in the reaction layout.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            layout (Reaction_Array): Well_Plate object containing the combinatorial reaction layout.

        Returns:
            dict: Dictionary containing the m/z values, retention times and smiles for the analytes.
        '''

        results_dict = {}
        for name, ms_injection in ms_sequence.injections.items():
            pos = ms_injection.get_plate_position()
            mz_match = ms_injection.match_mol(layout.get_product(pos), return_canidates=return_canidates, check_iso=check_iso, **kwargs)
            if mz_match:
                results_dict[pos] = [mz_match, layout.get_product(pos)]
            else:
                results_dict[pos] = [{}]
        return results_dict

