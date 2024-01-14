from pygecko import Agilent_MS_Parser
from pygecko.gc_tools import Spectral_Match
from pygecko.visualization import Visualization

# Load injection with pure analyte and pick peaks.
product_injection = Agilent_MS_Parser.load_injection('24a/FKB-FA-062-Prod.D')
product_injection.pick_peaks()

# Get the peak of the pure analyte.
product_peak = product_injection[9.409]

# Load sequence with plate data.
plate_sequence = Agilent_MS_Parser.load_sequence('Buchwald-Hartwig_Plate_GC_Data/MS', pos=True)

# Get the injection with the target product and pick peaks.
D1 = plate_sequence['FBS-FA-034-D1']
D1.pick_peaks()

# Find the peak in the injection that matches the pure analyte peak.
match = Spectral_Match.find_peak(product_peak, D1)

# Visualize the match.
Visualization.compare_mass_spectra(match)