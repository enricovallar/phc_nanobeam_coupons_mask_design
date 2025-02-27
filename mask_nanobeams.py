#________________________________________________________________
#       MASK FOR NANOBEMS WITH CAVITY
#________________________________________________________________
"""
    This script is used to draw a mask for uTP coupons
    with a nanobeam. 
    The nanobeam cavity parameters are swept as 
    well as the theter system parameters.
    The mask is exported to GDS.

    A first group has a central straight waveguide 30 um long, 
    a second group has a central straight waveguide 20 um long.
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
nanobeam_coupon = nanobeam_devices.NanobeamCoupon().build()
nanobeam_coupon.print_parameters()

nanobeam_coupon.set_param("theter_length", 1)
nanobeam_coupon.set_param("nanobeam_straight_length", 10)
nanobeam_coupon.set_param("nanobeam_taper_length", 8)
nanobeam_coupon.set_param("nanobeam_width", 1)
nanobeam_coupon.set_param("nanobeam_taper_m", 0.85)
nanobeam_coupon.set_param("theter_distance_from_tip", 8)
nanobeam_coupon.set_param("shrink_factor", 0.012)
nanobeam_coupon.set_param("r_c", 0.080)
nanobeam_coupon.set_param("n_left", 14)
nanobeam_coupon.set_param("theter_width", 0.2)
nanobeam_coupon.set_param("theter_space", 8)
nanobeam_coupon.set_param("theter_both_sides", 0)
nanobeam_coupon.set_param("shrink_um", 0.012)

## Sweep parameters
theter_space = [5, 10]
lattice_constant = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
lattice_constant = [0.350 + x for x in lattice_constant]
print("cavity lattice constants: ", lattice_constant)
gap = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
gap = [0.880 + x for x in gap]
print("cavity gaps: ", gap)
N_atoms_right = [14, 4]

param_cols, N_cols = generate_param_dict(
    a_c = lattice_constant, 
    n_right = N_atoms_right,
)

param_rows, N_rows = generate_param_dict(
    g_c = gap, 
)

row_id = [i for i in range(N_rows)] 
col_id = [i for i in range(N_cols)] 

param_rows["row_id"] = row_id
param_cols["col_id"] = col_id

print(param_cols)
print(param_rows)

## First group with 30um long nanobeam
nanobeam_coupon.set_param("nanobeam_straight_length", 30)   
nanobeam_coupon.set_param("group_id", 1)

# set annotations
nanobeam_coupon.annotate_params(
    "a_c",
    "n_right",
    "g_c",
    "nanobeam_straight_length",
)

print("Starting to draw the devices")
tab30um = smlay.DeviceTable(
    nanobeam_coupon, 
    N_rows,
    N_cols,
    param_rows,
    param_cols,
)
print("tab30um ready")

## Second group with 20um long nanobeam


nanobeam_coupon.set_param("nanobeam_straight_length", 20)
nanobeam_coupon.set_param("group_id", 2)

# set annotations
nanobeam_coupon.annotate_params(
    "a_c",
    "n_right",
    "g_c",
    "nanobeam_straight_length",
)

tab20um = smlay.DeviceTable(
    nanobeam_coupon,
    N_rows,
    N_cols,
    param_rows,
    param_cols,
)
print("tab20um ready")

tab30um.auto_align(20, 20, numkey=5)
tab20um.auto_align(20, 20, numkey=5)

print("Adding tables to the mask")
mask.addDeviceTable(tab30um, 0, 0)
mask.addDeviceTable(tab20um, 0, -300)
print("Exporting mask")
mask.exportGDS()






