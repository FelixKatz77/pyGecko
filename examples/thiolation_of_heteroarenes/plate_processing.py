from pygecko.parsers import Agilent_FID_Parser, Agilent_MS_Parser
from pygecko.gc_tools import load_sequence, save_sequence
from pygecko.reaction import Transformation, Well_Plate
from pygecko.visualization.visuals import Visualization
from pygecko.analysis.analysis import Analysis
from pygecko.reaction import Reaction_Parser

def main():

    # Files with GC-MS and GC-FID data of the reaction array.
    ms_path = 'Thiolation_Plate_GC_Data/MS'
    fid_path = 'Thiolation_Plate_GC_Data/FID'

    # Files and folders with the retention index calibration data.
    ms_ri_path = 'Thiolation_Plate_GC_Data/RI/MS'
    fid_ri_path = 'Thiolation_Plate_GC_Data/RI/FID'

    # Files with the reaction layout and reaction metadata.
    layout_path = 'Thiolation_Plate_GC_Data/plate_layout.csv'
    meta_data_path = 'Thiolation_Plate_GC_Data/meta_data.json'

    # Reaction SMARTS used to map substrates to products.
    rxn = Transformation(
        '[c:1]1([a:2][a:3][a:4][a:5]1)[Cl,Br:6].[SH1:7][#6:8]>>[c:1]1([a:2][a:3][a:4][a:5]1)[S:7][#6:8]')

    # Create a Well_Plate object from the reaction layout and reaction SMARTS.
    layout = Well_Plate(layout_path, rxn, meta_data_file=meta_data_path)

    # Load the GC-MS and GC-FID data.
    fid_sequence = load_sequence(fid_path)
    ms_sequence = load_sequence(ms_path)

    # Pick peaks in the GC-MS and GC-FID data.
    fid_sequence.pick_peaks()
    ms_sequence.pick_peaks()

    # Set the internal standard for the GC-FID data.
    fid_sequence.set_internal_standard(3.407, name='Dodecane', smiles='CCCCCCCCCCCC')

    # Load the retention index calibration data.
    ri_conf_ms = Agilent_MS_Parser.load_ri_calibration(ms_ri_path,12, rt=2.255)
    ri_conf_fid = Agilent_FID_Parser.load_ri_calibration(fid_ri_path, 2.4, c_count=12, rt=3.394)

    # Assign retention indices to the GC-MS and GC-FID data.
    ri_conf_ms.assign_ris(ms_sequence)
    ri_conf_fid.assign_ris(fid_sequence, alignment=True)

    # Calculate the yields of the reactions.
    yield_array = Analysis.calc_plate_yield(ms_sequence, fid_sequence, layout)

    # Generate plate heatmap.
    Visualization.visualize_plate(yield_array['quantity'], well_labels=True,
                                  row_labels=["1", "2", "3", "4", "5", "6$^a$", "7$^a$", "8"],
                                  col_labels=["9", "10", "11", "12", "13", "14", "15", "16", "17", "18$^b$", "19",
                                              "20"])
    # Generate ORD dataset.
    Reaction_Parser.build_dataset(layout, yield_array['quantity'], 'heteroarene_thiolation.pbtxt')

if __name__ == '__main__':
    main()