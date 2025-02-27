#________________________________________________________________
##      MASK FOR NANOBEMS WITH PHC
#________________________________________________________________
""" This script is used to draw a mask for uTP coupons
with a nanobeam realized inside a PhC.

First we sweep the parameters for the nanobeam cavity, 
then for the theter system, without including the cavity.
The mask is exported to GDS.
"""

# Let's import basic stuff
import samplemaker.layout as smlay # used for layout 
import samplemaker.makers as sm # used for drawing
import samplemaker.devices as smdev # used for device function
# Let's use numpy arrays
import numpy as np

# Let's import the custom devices
import nanobeam_devices
from nanobeam_devices import *


# Get the mask name from the file name of this script
import os
mask_name = os.path.splitext(os.path.basename(__file__))[0]
mask = smlay.Mask(mask_name)
mask.set_cache(True)
mask.addWriteFieldGrid(500, 0, 0, 2, 2)


# e-beam markers
markdev = smdev.Device.build_registered("BASELIB_CMARK")
markerset = smlay.MarkerSet("Ebeam1", markdev, x0=-200, y0=-200, mset=4, xdist=900, ydist=900)
mask.addMarkers(markerset)

# Add nano-beam devices
phc_nb_coupon = nanobeam_devices.PhCNanobeamCoupon().build()

#DeviceInspect(phc_nb_coupon) # ------UNCOMMENT TO INSPECT DEVICE
phc_nb_coupon.print_parameters()

# set fixed parameters
phc_nb_coupon.set_param("theter_length", 1)
phc_nb_coupon.set_param("nanobeam_taper_length", 8)
phc_nb_coupon.set_param("nanobeam_width", 1)
phc_nb_coupon.set_param("nanobeam_taper_m", 0.85)
phc_nb_coupon.set_param("theter_distance_from_corner", 2)
phc_nb_coupon.set_param("shrink_um", 0.012)
phc_nb_coupon.set_param("r_c", 0.080)
phc_nb_coupon.set_param("n_left", 14)
phc_nb_coupon.set_param("theter_width", 0.2)
phc_nb_coupon.set_param("theter_separation", 8)
phc_nb_coupon.set_param("pad_edge_x", 40)
phc_nb_coupon.set_param("pad_edge_y", 10)
phc_nb_coupon.set_param("nanobeam_straight_length", 20)


# build phc referenced later (required for increased writing speed)
phc = nanobeam_devices.make_phc_reference(length = 35, Ny=5, 
                                          lattice_constant=0.455, radius=0.11,
                                          mask=mask)
phc_nb_coupon.set_mirror_phc(phc)


#---------------------------------------------------------
## Sweep parameters for the cavity
#---------------------------------------------------------
lattice_constant = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
lattice_constant = [0.350 + x for x in lattice_constant]
print("cavity lattice constants: ", lattice_constant)
gap = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
gap = [0.880 + x for x in gap]
print("cavity gaps: ", gap)
N_atoms_right = [14,  4]
N_atoms_left  = [14]

param_cols, N_cols = generate_param_dict(
    a_c = lattice_constant, 
    n_right = N_atoms_right,
)

param_rows, N_rows = generate_param_dict(
    g_c = gap, 
    n_left = N_atoms_left,
)

phc_nb_coupon.annotate_params("a_c",
                              "g_c",
                              "n_right",
                              "n_left",)

row_id = [i for i in range(N_rows)] 
col_id = [i for i in range(N_cols)] 

param_rows["row_id"] = row_id
param_cols["col_id"] = col_id


phc_nb_coupon.set_param("nanobeam_straight_length", 20)   
phc_nb_coupon.set_param("group_id", 3)

print("Starting to draw the devices")
tab_cavity = smlay.DeviceTable(
    phc_nb_coupon, 
    N_rows,
    N_cols,
    param_rows,
    param_cols,
)
print("tab_cavity ready")

#---------------------------------------------------------
## Sweep parameters for the thetes, no cavity
#---------------------------------------------------------

theter_width = np.arange(0.150, 0.35, 0.05).tolist()    
theter_separation = [4, 6, 8]
nanobeam_straight_length = [20, 20] 

param_cols, N_cols = generate_param_dict(
    theter_width = theter_width,
    nanobeam_straight_length = nanobeam_straight_length 
)

param_rows, N_rows = generate_param_dict(
    theter_separation = theter_separation,
)

phc_nb_coupon.annotate_params(
    "theter_width",
    "nanobeam_straight_length",
    "theter_separation",
)

row_id = [i for i in range(N_rows)] 
col_id = [i for i in range(N_cols)] 

param_rows["row_id"] = row_id
param_cols["col_id"] = col_id

phc_nb_coupon.set_param("n_left", 0)
phc_nb_coupon.set_param("n_right", 0)
phc_nb_coupon.set_param("group_id", 4)

tab_no_cavity = smlay.DeviceTable(
    phc_nb_coupon,
    N_rows,
    N_cols,
    param_rows,
    param_cols,
)
print("tab_no_cavity ready")

## Write the mask  to GDS
# Align the tables
tab_cavity.auto_align(20, 20, numkey=5)
tab_no_cavity.auto_align(20, 20, numkey=5)

# Add tables to the mask
print("Adding tables to the mask")
mask.addDeviceTable(tab_cavity, 0, 0)
mask.addDeviceTable(tab_no_cavity, 0, -300)
print("Exporting mask")
mask.exportGDS()






