import numpy as np
import pandas as pd
from pygecko.gc_tools.sequence import MS_Sequence, FID_Sequence
from pygecko.gc_tools.analyte import Analyte
from pygecko.reaction import Well_Plate
from numpy.lib.recfunctions import unstructured_to_structured


class Analysis:

    '''
    A class wrapping functions to analyze GC data.
    '''

    @staticmethod
    def calc_plate_yield(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Well_Plate,
                                 path: str|None = None):

        '''
        Matches GC-MS and GC-FID peaks and quantifies the yields of the reactions.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.
            path (str|None, optional): Path to write the results to. Defaults to None.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''

        return Analysis.__match_and_quantify_plate(ms_sequence, fid_sequence, layout, path, mode='yield')

    @staticmethod
    def calc_plate_conv(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Well_Plate,
                                 path: str|None = None, index:int=0):

        '''
        Matches GC-MS and GC-FID peaks and quantifies the conversion of the reactions.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.
            path (str|None, optional): Path to write the results to. Defaults to None.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''

        return Analysis.__match_and_quantify_plate(ms_sequence, fid_sequence, layout, path, mode='conv', index=index)

    @staticmethod
    def __match_and_quantify_plate(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Well_Plate,
                                 path: str|None = None, mode:str='yield', index:int=0) -> np.ndarray:

        '''
        Matches GC-MS and GC-FID peaks and quantifies the yields/conversion of the reactions.

        Args:
            ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.
            path (str|None, optional): Path to write the results to. Defaults to None.
            mode (str, optional): Parameter (yield or conversion) to quantify. Defaults to 'yield'.

        Returns:
            np.ndarray: Numpy array containing the quantification results, retention times and smiles for the analytes.
        '''

        results_df = pd.DataFrame(columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
                                      index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

        results_dict = Analysis.__match_and_quantify(ms_sequence, fid_sequence, layout, mode, index)

        for key, value in results_dict.items():
            results_df.loc[key[0], key[1:]] = value

        results_array = results_df.to_numpy()
        results_array = np.array([row.tolist() for row in results_array])
        dtype = np.dtype([('quantity', float), ('rt_ms', float), ('rt_fid', float)])
        results_array = unstructured_to_structured(results_array[:,:,:-1], dtype=dtype)
        if path:
            if mode == 'yield':
                quantity = 'Yield [%]'
                report_df = pd.DataFrame.from_dict(results_dict, orient='index', columns=[quantity, 'RT-MS [min]', 'RT-FID [min]', 'Analyte'])
            else:
                quantity = 'Conversion [%]'
                report_df = pd.DataFrame.from_dict(results_dict, orient='index', columns=[quantity, 'RT-MS [min]', 'RT-FID [min]', 'Analyte'])
            report_df.sort_index(inplace=True)
            s = report_df.style.background_gradient(axis=0, subset=quantity, cmap='GnBu')
            s.applymap(lambda x: 'color: transparent; background-color: transparent' if pd.isnull(x) else '').to_excel(path, engine='openpyxl')
        return results_array

    @staticmethod
    def __match_and_quantify(ms_sequence: MS_Sequence, fid_sequence: FID_Sequence, layout: Well_Plate, mode: str,
                           index: int = 0) -> dict[str, list[float, str]]:

        '''
                Matches GC-MS and GC-FID peaks and quantifies the yields of the reactions.

                Args:
                    ms_sequence (MS_Sequence): MS_Sequence object containing the GC-MS data.
                    fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
                    layout (Well_Plate): Well_Plate object containing the combinatorial reaction layout.

                Returns:
                    dict: Dictionary containing the yields, retention times and the analyte smiles for each well.
        '''

        results_dict = {}
        for name, ms_injection in ms_sequence.injections.items():
            fid_injection = fid_sequence[name]
            pos = ms_injection.get_plate_position()
            if pos == 'F3':
                pass
            if mode == 'yield':
                analyte = layout.get_product(ms_injection.get_plate_position())
                pass
            if mode == 'conv':
                analyte = layout.get_substrate(pos, index=index)
            mz_match = ms_injection.match_mol(analyte)

            if mz_match:
                ri_match = fid_injection.match_ri(mz_match.ri, analyte=mz_match.analyte)
                if ri_match:
                    yield_ = fid_injection.quantify(ri_match.rt)
                    if mode == 'conv':
                        yield_ = 100 - yield_
                    results_dict[pos] = [yield_, mz_match.rt, ri_match.rt, analyte]
                else:
                    results_dict[pos] = [np.nan, np.nan, np.nan, '']
            else:
                results_dict[pos] = [np.nan, np.nan, np.nan, '']
        return results_dict


    @staticmethod
    def quantify_analyte(fid_sequence:FID_Sequence, rt:float, analyte:Analyte|None=None) -> dict[str, float]:

        '''
        Returns the yields of Analytes at a given retention time in a FID sequence.

        Args:
            fid_sequence (FID_Sequence): FID_Sequence object containing the GC-FID data.
            rt (float): Retention time to quantify the Analytes at.
            analyte (Analyte): Analyte to quantify.

        Returns:
            dict: Dictionary containing the yields of the Analytes in the FID sequence.
        '''

        yield_dict = {}
        for injection in fid_sequence.injections.values():
            peak = injection.flag_peak(rt, analyte=analyte)
            if peak:
                yield_ = injection.quantify(peak.rt)
            else:
                yield_ = 0
            if injection.plate_pos:
                yield_dict[injection.plate_pos] = yield_
            else:
                yield_dict[injection.sample_name] = yield_
        return yield_dict