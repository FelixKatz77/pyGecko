from pygecko.parsers import Agilent_FID_Parser, MS_Base_Parser
from pygecko.reaction import Transformation, Reaction_Array
from pygecko.visualization.visuals import Visualization
from pygecko.analysis.analysis import Analysis
from pygecko.reaction import Reaction_Parser

def main():

    # Files with GC-MS and GC-FID data of the reaction array.
    ms_path = 'Buchwald-Hartwig_Plate_GC_Data/MS'
    fid_path = 'Buchwald-Hartwig_Plate_GC_Data/FID'

    # Files and folders with the retention index calibration data.
    ms_ri_path = 'Buchwald-Hartwig_Plate_GC_Data/RI/MS/FBS-FA-033-RI-II.mzML'
    fid_ri_path = 'Buchwald-Hartwig_Plate_GC_Data/RI/FID'

    # Files with the reaction layout and reaction metadata.
    layout_path = 'Buchwald-Hartwig_Plate_GC_Data/plate_layout.csv'
    meta_data_path = 'Buchwald-Hartwig_Plate_GC_Data/meta_data.json'

    # Reaction SMARTS used to map substrates to products.
    rxn = Transformation(
        '[C,c:1][Nh1,Nh2,nh1:2].[Br,Cl:3][C,c:4]>>[C,c:1][N,n:2][C,c:4]')

    # Create a Well_Plate object from the reaction layout and reaction SMARTS.
    layout = Reaction_Array(layout_path, rxn, meta_data_file=meta_data_path)

    # Load the GC-MS and GC-FID data.
    fid_sequence = Agilent_FID_Parser.load_sequence(fid_path, 2.7, pos=True)
    ms_sequence = MS_Base_Parser.load_sequence(ms_path, pos=True)

    # Pick peaks in the GC-MS and GC-FID data.
    fid_sequence.pick_peaks()
    ms_sequence.pick_peaks(prominence_ms=125)

    # Set the internal standard.
    fid_sequence.set_internal_standard(4.593, name='Dodecane', smiles='CCCCCCCCCCCC')
    ms_sequence.set_internal_standard(3.324, name='Dodecane', smiles='CCCCCCCCCCCC')

    # Load the retention index calibration data.
    ri_conf_ms = MS_Base_Parser.load_ri_calibration(ms_ri_path, 10, rt=2.154)
    ri_conf_fid = Agilent_FID_Parser.load_ri_calibration(fid_ri_path, 2.7, c_count=12, rt=4.593)

    # Assign retention indices to the GC-MS and GC-FID data.
    ri_conf_ms.assign_ris(ms_sequence)
    ri_conf_fid.assign_ris(fid_sequence, alignment=True)

    # Calculate the yields of the reactions.
    yield_array = Analysis.calc_plate_yield(ms_sequence, fid_sequence, layout)

    # Generate plate heatmap.
    Visualization.visualize_plate(yield_array, well_labels=True,
                                  row_labels=["21", "22", "23", "24", "25", "26", "27", "28"],
                                  col_labels=["29$^a$", "30$^a$", "31$^a$", "32$^a$", "33$^a$", "34$^a$", "29$^b$",
                                              "30$^b$", "31$^b$", "32$^b$", "33$^b$", "34$^b$"])
    # Generate ORD dataset.
    Reaction_Parser.build_dataset(layout, yield_array['quantity'], 'buchwald_hartwig.pbtxt')

if __name__ == '__main__':
    main()