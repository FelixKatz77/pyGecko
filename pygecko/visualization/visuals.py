import numpy as np
from typing import Union
from typing import TYPE_CHECKING
from pygecko.gc_tools.utilities import Utilities
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Wedge
import matplotlib
from matplotlib.ticker import (MultipleLocator)
from matplotlib.figure import Figure, figaspect
from pygecko.visualization.utilities import yield_cmap

if TYPE_CHECKING:
    from pygecko.gc_tools import Injection, FID_Injection, MS_Injection, MS_Peak

# Applied per figure through plt.rc_context, never to the global rcParams: setting those at import
# restyled every other matplotlib user in the process.
_STYLE = {'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 12,
          'font.weight': 'regular'}

class Visualization:

    @staticmethod
    def visualize_plate(data: np.ndarray, well_labels=True, show_flags:bool=False, **kwargs) -> Figure:

        '''
        Returns a well plate visualized as a heatmap of yields.

        Args:
            data (np.ndarray): A numpy array containing the yields of the reactions.
            cbar_label (str, optional): Label for the colorbar. Defaults to 'Yield [%]'; pass e.g.
                'Conversion [%]' when plotting conversion data.
        '''

        row_labels = kwargs.pop('row_labels', ["A", "B", "C", "D", "E", "F", "G", "H"])
        col_labels = kwargs.pop('col_labels', ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"])
        cbar_label = kwargs.pop('cbar_label', 'Yield [%]')

        flags = data['flags']
        data = data['quantity']
        masked_data = np.ma.array (data, mask=np.isnan(data))
        cmap, norm = yield_cmap
        cmap.set_bad('darkgrey', 0.5)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
        N = data.shape[0]
        M = data.shape[1]
        r = 0.43

        x, y = np.meshgrid(np.arange(M), np.arange(N))

        with plt.rc_context(_STYLE):
            fig, ax = plt.subplots(figsize=(8.5, 4.8))
            ax.invert_yaxis()
            circles = [Circle((j, i), radius=r) for j, i in zip(x.flat, y.flat)]
            col = PatchCollection(circles, array=masked_data.flatten(), cmap=cmap, norm=norm)
            ax.add_collection(col)

            if show_flags:
                for i in range(N):
                    for j in range(M):
                        if flags[i, j] == 1:
                            indicator = Wedge((j+0.3, i-0.3), r * 0.4, 0, 360, color='orange')
                            ax.add_patch(indicator)
                            ax.text(j + 0.3, i - 0.29, '!', color='#a44018', fontsize=12,
                                    ha='center', va='center', fontweight='bold')

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
            cbar = fig.colorbar(sm, ax=ax, ticks=[0, 25, 50, 75, 100])
            cbar.ax.set_ylabel(cbar_label, size=14)
            cbar.ax.tick_params(labelsize=12, )
            fig.tight_layout()
        return fig



    @staticmethod
    def view_chromatogram(injection:Union['MS_Injection', 'FID_Injection'], **kwargs) -> Figure:

        '''
        Returns a chromatogram visualized as time/intensity plot.

        Args:
            injection ('MS_Injection'|'FID_Injection'): Injection object containing the chromatogram to visualize.
            **kwargs: Keyword arguments for the plot.
        '''

        raw = kwargs.pop('raw', False)
        chromatogram = injection.chromatogram
        x = chromatogram.time
        y = chromatogram.intensity if raw or injection.detector == 'MS' else chromatogram.processed

        xlim = kwargs.pop('xlim', [x.min(), x.max()])
        ylim = kwargs.pop('ylim', [-100000, None])
        highlight_peaks = kwargs.pop('highlight_peaks', True)
        color = kwargs.pop('color', '#005573')
        linewidth = kwargs.pop('linewidth', 1)

        w, h = figaspect(0.5)
        with plt.rc_context(_STYLE):
            fig, ax = plt.subplots(figsize=(w, h))
            xlim_scans = Utilities.convert_time_to_scan([x - injection.solvent_delay for x in xlim],
                                                        chromatogram.scan_rate)

            x, y = x[xlim_scans[0]:xlim_scans[1]], y[xlim_scans[0]:xlim_scans[1]]

            ax.plot(x, y, color=color, lw=linewidth, **kwargs)

            if highlight_peaks:
                ax.fill_between(x, y.min(), y.max(), where=injection._check_for_peak(x),
                                 color='#005573', alpha=0.1, transform=ax.get_xaxis_transform())


            ax.grid(color='lightgrey', linestyle='--', which='both')

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.25))

            ax.set_xlabel('Time [min]')
            ax.set_title(injection.sample_name)

            fig.tight_layout()
        return fig

    @staticmethod
    def view_mass_spectrum(peak:'MS_Peak', **kwargs) -> Figure:

        '''
        Returns a mass spectrum visualized as m/z/intensity plot.

        Args:
            peak (MS_Peak): MS_Peak object containing the mass spectrum to visualize.
            **kwargs: Keyword arguments for the plot.
        '''

        xlim = kwargs.pop('xlim', [None, None])
        ylim = kwargs.pop('ylim', [None, None])
        color = kwargs.pop('color', '#005573')

        with plt.rc_context(_STYLE):
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

            ax.grid(color='lightgrey', linestyle='--', which='both', axis='y')
            ax.set_xlabel('m/z')
            ax.set_ylabel('Relative Intensity [%]')
            ax.spines[['right', 'top']].set_visible(False)
            fig.tight_layout()

        return fig

    @staticmethod
    def stack_chromatograms(injections:list['Injection'], **kwargs) -> Figure:

        '''
        Returns a list of chromatograms visualized as a stack plot.

        Args:
            injections (list['Injection']): List of Injection objects containing the chromatograms to visualize.
            **kwargs: Keyword arguments for the plot.
        '''

        color = kwargs.pop('color', '#005573')
        linewidth = kwargs.pop('linewidth', 1)
        # Popped once, before the loop, so they apply to every trace rather than to the first only.
        raw = kwargs.pop('raw', False)
        xlim = kwargs.pop('xlim', None)
        ylim = kwargs.pop('ylim', [-100000, None])

        with plt.rc_context(_STYLE):
            fig, axes = plt.subplots(len(injections), 1, figsize=(8, 6), squeeze=False)

            ax_objs = []

            for i, injection in enumerate(injections):

                if isinstance(color, list) and len(color) == len(injections):
                    color_i = color[i]
                elif isinstance(color, str):
                    color_i = color
                else:
                    color_i = '#005573'

                ax = axes[i, 0]

                chromatogram = injection.chromatogram
                x = chromatogram.time
                y = chromatogram.intensity if raw or injection.detector == 'MS' else chromatogram.processed

                xlim_i = xlim if xlim is not None else [x.min(), x.max()]

                xlim_scans = Utilities.convert_time_to_scan([x - injection.solvent_delay for x in xlim_i],
                                                            chromatogram.scan_rate)
                x, y = x[xlim_scans[0]:xlim_scans[1]], y[xlim_scans[0]:xlim_scans[1]]

                ax.plot(x, y, color=color_i, lw=linewidth, **kwargs)
                ax.set_xlim(xlim_i)
                ax.set_ylim(ylim)

                ax.xaxis.set_major_locator(MultipleLocator(1.0))
                ax.xaxis.set_minor_locator(MultipleLocator(0.5))
                ax.grid(color='lightgrey', linestyle='--', which='both', axis='x')
                if i == len(injections) - 1:
                    ax.set_xlabel('Time [min]')
                if i != len(injections) - 1:
                    ax.set_xticklabels([])
                ax.set_ylabel(f'{injection.sample_name}')
                ax.spines[['right', 'top']].set_visible(False)

                ax_objs.append(ax)

            fig.tight_layout()

        return fig

    @staticmethod
    def compare_mass_spectra(peaks:tuple['MS_Peak', 'MS_Peak'], **kwargs) -> Figure:

        '''
        Returns two mass spectra visualized as m/z/intensity plot.

        Args:
            peaks(tuple['MS_Peak', 'MS_Peak']): Tuple of MS_Peak objects containing the mass spectra to visualize.
            **kwargs: Keyword arguments for the plot.
        '''

        peak1, peak2 = peaks
        xlim = kwargs.pop('xlim', [None, None])
        ylim = kwargs.pop('ylim', [None, None])
        colors = kwargs.pop('colors', ('#005573', '#e04214'))

        with plt.rc_context(_STYLE):
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

            ax.grid(color='lightgrey', linestyle='--', which='both', axis='y')
            ax.set_xlabel('m/z')
            ax.set_ylabel('Relative Intensity [%]')
            ax.set_yticks([-100, -75, -50, -25, 0, 25, 50, 75, 100])
            ax.set_yticklabels([100, 75, 50, 25, 0, 25, 50, 75, 100])
            ax.spines[['right', 'top']].set_visible(False)
            fig.tight_layout()

        return fig

