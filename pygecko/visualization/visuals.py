import numpy as np
from typing import Union
from typing import TYPE_CHECKING
from pygecko.gc_tools.utilities import Utilities
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib
from matplotlib.ticker import (MultipleLocator)
from matplotlib.figure import figaspect
import matplotlib.gridspec as grid_spec
from pygecko.visualization.utilities import yield_cmap

if TYPE_CHECKING:
    from pygecko.gc_tools import Injection, FID_Injection, MS_Injection, MS_Peak

plt.rcParams["font.family"] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 12
plt.rcParams['font.weight'] = 'regular'

class Visualization:

    @staticmethod
    def visualize_plate(data: np.ndarray, path:str|None=None, well_labels=True, **kwargs) -> None:

        '''
        Visualizes a well plate as a heatmap of yields and saves the figure if a path is given.

        Args:
            data (np.ndarray): A numpy array containing the yields of the reactions.
            results (str, optional): The type results to visualize. Defaults to 'hit'.
            path (str|None, optional): Path to save the figure to. Defaults to None.
        '''

        row_labels = kwargs.pop('row_labels', ["A", "B", "C", "D", "E", "F", "G", "H"])
        col_labels = kwargs.pop('col_labels', ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"])

        masked_data = np.ma.array (data, mask=np.isnan(data))
        cmap, norm = yield_cmap
        cmap.set_bad('darkgrey', 0.5)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
        N = data.shape[0]
        M = data.shape[1]
        r = 0.43

        x, y = np.meshgrid(np.arange(M), np.arange(N))

        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        plt.gca().invert_yaxis()
        circles = [plt.Circle((j, i), radius=r) for j, i in zip(x.flat, y.flat)]
        col = PatchCollection(circles, array=masked_data.flatten(), cmap=cmap, norm=norm)
        ax.add_collection(col)

        ax.set_xticks(np.arange(data.shape[1]), labels=col_labels, weight='bold')
        ax.set_yticks(np.arange(data.shape[0]), labels=row_labels, weight='bold')
        ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
        ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
        ax.tick_params(top=True, bottom=False,
                       labeltop=True, labelbottom=False, length=0)
        ax.tick_params(axis='y', which='major', pad=7)
        ax.spines[:].set_visible(False)

        ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
        ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
        ax.grid(which="minor", color="darkgrey", linestyle='-', linewidth=1)
        ax.tick_params(which="minor", bottom=False, left=False)

        fig.patch.set_facecolor('white')
        fig.patch.set_alpha(0.0)
        ax.set_facecolor('lightgrey')

        if well_labels:
            for i in range(N):
                for j in range(M):
                    if not np.isnan(data[i, j]):
                        ax.text(j, i, int(round(data[i, j], 0)), ha="center", va="center", color="black")


        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ticks=[0, 25, 50, 75, 100])
        cbar.ax.set_ylabel('Yield [%]', size=14)
        cbar.ax.tick_params(labelsize=12, )
        fig.tight_layout()
        if path:
            plt.savefig(path, dpi=400)
            plt.close()
        else:
            plt.show()



    @staticmethod
    def view_chromatogram(injection:Union['MS_Injection', 'FID_Injection'], path:str|None=None, **kwargs):

        '''
        Visualizes a chromatogram as time/intensity plot.

        Args:
            injection ('MS_Injection'|'FID_Injection'): Injection object containing the chromatogram to visualize.
            path (str|None, optional): Path to save the figure to. Defaults to None.
            **kwargs: Keyword arguments for the plot.
        '''

        raw = kwargs.pop('raw', False)
        if raw or injection.detector == 'MS':
            chromatogram = injection.chromatogram
        else:
            chromatogram = injection.processed_chromatogram
        x = chromatogram[0]
        y = chromatogram[1]

        xlim = kwargs.pop('xlim', [x.min(), x.max()])
        ylim = kwargs.pop('ylim', [-100000, None])
        highlight_peaks = kwargs.pop('highlight_peaks', True)
        color = kwargs.pop('color', '#005573')
        linewidth = kwargs.pop('linewidth', 1)

        w, h = figaspect(0.5)
        fig, ax = plt.subplots(figsize=(w, h))
        xlim_scans = Utilities.convert_time_to_scan([x - injection.solvent_delay for x in xlim],
                                                    injection.analysis_settings.scan_rate)

        x, y = x[xlim_scans[0]:xlim_scans[1]], y[xlim_scans[0]:xlim_scans[1]]

        ax.plot(x, y, color=color, lw=linewidth, **kwargs)

        if highlight_peaks and injection.detector == 'FID':
            ax.fill_between(x, y.min(), y.max(), where=injection._check_for_peak(x),
                             color='#005573', alpha=0.1, transform=ax.get_xaxis_transform())


        plt.grid(color='lightgrey', linestyle='--', which='both')

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))

        ax.set_xlabel('Time [min]')
        plt.title(injection.sample_name)

        fig.tight_layout()
        if path:
            plt.savefig(path, dpi=400)
            plt.close()
        else:
            plt.show()

    @staticmethod
    def view_mass_spectrum(peak:'MS_Peak', path:str|None=None, **kwargs) -> None:

        '''
        Visualizes a mass spectrum as m/z/intensity plot.

        Args:
            peak (MS_Peak): MS_Peak object containing the mass spectrum to visualize.
            path (str|None, optional): Path to save the figure to. Defaults to None.
            **kwargs: Keyword arguments for the plot.
        '''

        xlim = kwargs.pop('xlim', [None, None])
        ylim = kwargs.pop('ylim', [None, None])
        color = kwargs.pop('color', '#005573')

        fig, ax = plt.subplots(figsize=(14.4, 4.8))
        ax.bar(peak.mass_spectrum['mz'], peak.mass_spectrum['rel_intensity'], width=0.05, color=color, edgecolor=color,
               **kwargs)

        lable_indices = np.argpartition(peak.mass_spectrum['rel_intensity'], -4)[-4:]
        for index in lable_indices:
            ax.annotate(f'{peak.mass_spectrum["mz"][index]:.0f}',
                        (peak.mass_spectrum['mz'][index], peak.mass_spectrum['rel_intensity'][index]),
                        textcoords="offset points", xytext=(0, 5), ha='center', color='darkgrey', fontsize=10)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        plt.grid(color='lightgrey', linestyle='--', which='both', axis='y')
        ax.set_xlabel('m/z')
        ax.set_ylabel('Relative Intensity [%]')
        ax.spines[['right', 'top']].set_visible(False)
        fig.tight_layout()

        if path:
            plt.savefig(path, dpi=400)
            plt.close()
        else:
            plt.show()

    @staticmethod
    def stack_chromatograms(injections:list['Injection'], path:str|None=None, **kwargs) -> None:

        '''
        Visualizes a list of chromatograms as a stack plot.

        Args:
            injections (list['Injection']): List of Injection objects containing the chromatograms to visualize.
            path (str|None, optional): Path to save the figure to. Defaults to None.
            **kwargs: Keyword arguments for the plot.
        '''

        color = kwargs.pop('color', '#005573')
        linewidth = kwargs.pop('linewidth', 1)

        gs = (grid_spec.GridSpec(len(injections), 1))
        fig = plt.figure(figsize=(8, 6))

        ax_objs = []

        for i, injection in enumerate(injections):

            if isinstance(color, list) and len(color) == len(injections):
                color_i = color[i]
            elif isinstance(color, str):
                color_i = color
            else:
                color_i = '#005573'

            ax = fig.add_subplot(gs[i:i+1, 0:])

            raw = kwargs.pop('raw', False)

            if raw or injection.detector == 'MS':
                chromatogram = injection.chromatogram
            else:
                chromatogram = injection.processed_chromatogram
            x = chromatogram[0]
            y = chromatogram[1]

            xlim = kwargs.pop('xlim', [x.min(), x.max()])
            ylim = kwargs.pop('ylim', [-100000, None])

            xlim_scans = Utilities.convert_time_to_scan([x - injection.solvent_delay for x in xlim],
                                                        injection.analysis_settings.scan_rate)
            x, y = x[xlim_scans[0]:xlim_scans[1]], y[xlim_scans[0]:xlim_scans[1]]

            ax.plot(x, y, color=color_i, lw=linewidth, **kwargs)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.xaxis.set_major_locator(MultipleLocator(1.0))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            plt.grid(color='lightgrey', linestyle='--', which='both', axis='x')
            if i == len(injections) - 1:
                ax.set_xlabel('Time [min]')
            if i != len(injections) - 1:
                ax.set_xticklabels([])
            ax.set_ylabel(f'{injection.sample_name}')
            ax.spines[['right', 'top']].set_visible(False)

            ax_objs.append(ax)

        plt.tight_layout()

        if path:
            plt.savefig(path, dpi=400)
            plt.close()
        else:
            plt.show()

    @staticmethod
    def compare_mass_spectra(peaks:tuple['MS_Peak', 'MS_Peak'], path:str|None=None, **kwargs) -> None:

        '''
        Visualizes two mass spectra as m/z/intensity plot.

        Args:
            peaks(tuple['MS_Peak', 'MS_Peak']): Tuple of MS_Peak objects containing the mass spectra to visualize.
            path (str|None, optional): Path to save the figure to. Defaults to None.
            **kwargs: Keyword arguments for the plot.
        '''

        peak1, peak2 = peaks
        xlim = kwargs.pop('xlim', [None, None])
        ylim = kwargs.pop('ylim', [None, None])
        colors = kwargs.pop('colors', ('#005573', '#e04214'))

        fig, ax = plt.subplots(figsize=(14.4, 4.8))
        ax.bar(peak1.mass_spectrum['mz'], peak1.mass_spectrum['rel_intensity'], width=0.05, color=colors[0],
               edgecolor=colors[0],
               **kwargs)
        ax.bar(peak2.mass_spectrum['mz'], peak2.mass_spectrum['rel_intensity']*(-1), width=0.05, color=colors[1],
               edgecolor=colors[1],
               **kwargs)

        lable_indices1 = np.argpartition(peak1.mass_spectrum['rel_intensity'], -4)[-4:]
        lable_indices2 = np.argpartition(peak2.mass_spectrum['rel_intensity'], -4)[-4:]

        for index in lable_indices1:
            ax.annotate(f'{peak1.mass_spectrum["mz"][index]:.0f}',
                        (peak1.mass_spectrum['mz'][index], peak1.mass_spectrum['rel_intensity'][index]),
                        textcoords="offset points", xytext=(0, 5), ha='center', color='darkgrey', fontsize=10)

        for index in lable_indices2:
            ax.annotate(f'{peak2.mass_spectrum["mz"][index]:.0f}',
                        (peak2.mass_spectrum['mz'][index], peak2.mass_spectrum['rel_intensity'][index]*(-1)),
                        textcoords="offset points", xytext=(0, -10), ha='center', color='darkgrey', fontsize=10)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        plt.grid(color='lightgrey', linestyle='--', which='both', axis='y')
        ax.set_xlabel('m/z')
        ax.set_ylabel('Relative Intensity [%]')
        ax.set_yticks([-100, -75, -50, -25, 0, 25, 50, 75, 100])
        ax.set_yticklabels([100, 75, 50, 25, 0, 25, 50, 75, 100])
        ax.spines[['right', 'top']].set_visible(False)
        fig.tight_layout()

        if path:
            plt.savefig(path, dpi=400)
            plt.close()
        else:
            plt.show()


if __name__ == '__main__':
    array_8x12 = np.random.randint(0, 95, size=(8, 12))
    Visualization.visualize_plate(array_8x12, path='mock_heatmap.png')
    pass