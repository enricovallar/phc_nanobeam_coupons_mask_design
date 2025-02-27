import samplemaker.layout as smlay
import samplemaker.makers as sm
from samplemaker.devices import Device, registerDevicesInModule
from samplemaker.viewers import DeviceInspect, GeomView
from samplemaker.baselib.waveguides import BaseWaveguideSequencer, BaseWaveguidePort
from samplemaker.routers import WaveguideConnect
from samplemaker.shapes import GeomGroup
# Create a simple mask layout
import numpy as np
import math
import samplemaker.devices


# Class cavity: a simple cavity 

class DeviceWithNotes(Device):  
    def __init__(self):
        super().__init__()
        self.relevant_params = []

    def initialize(self):
        self.set_name("DEVICE_WITH_NOTES")
        self.set_description("A simple device that shows the parameters in a text box")   
    
    def write_notes(self):
        p = self.get_params()
        t = ""
        for key in self.relevant_params:
            t += f"{key}: {p[key]: 0.3f}\n"
        text = sm.make_text(0,0, t, 1, 0.2, to_poly=False, numkey=8, layer=7)
        return text
    
    def print_parameters(self):
        p = self.get_params()
        print(f"{self.name} parameters and default values:")
        for key, value in p.items():
            print(f"  {key}: {value}")

    def parameters(self):
        self.addparameter("row_id", 1, "Row ID", param_type=int)
        self.addparameter("col_id", 1, "Column ID", param_type=int)
        self.addparameter("group_id", 1, "Group ID", param_type=int)

    def geom(self):
        p = self.get_params()
        
        id = f"{p['group_id']}_{p['row_id']}_{p['col_id']}"
        id_text = sm.make_text(0, 0, id, 2, 0.2, to_poly=True, numkey=2, layer=9)
        return id_text
    
    def make_id_text(self):
        p = self.get_params()
        id = f"{p['group_id']}_{p['row_id']}_{p['col_id']}"
        id_text = sm.make_text(0, 0, id, 2, 0.2, to_poly=True, numkey=2, layer=9)
        return id_text



    @property
    def name(self):
        return self._name
    
    @property
    def relevant_params(self):
        return self._relevant_params
    
    @relevant_params.setter
    def relevant_params(self, relevant_params):
        self._relevant_params = relevant_params

    def annotate_params(self, *params):
        self.relevant_params = params

    


class Cavity(DeviceWithNotes):
    def __init__(self):
        super().__init__()
        self.relevant_params = []

    def initialize(self):
        self.set_name("CAVITY")
        self.set_description("A simple photonic crystal cavity")

    def parameters(self):
        super().parameters()
        self.addparameter("r_c", 0.080, "Radius of the holes of the cavity")
        self.addparameter("n_left", 5, "Number of holes on the left side", param_type=int)
        self.addparameter("n_right", 5, "Number of holes on the right side", param_type=int)
        self.addparameter("g_c" , 0.880, "Gap between the centers of the two holes of the cavity")
        self.addparameter("a_c", 0.350, "Lattice constant of the cavity")
    
    
    def geom(self): 
        p = self.get_params()
        x_centers_left = []
        x = 0
        for i in range(p["n_left"]):
            x_centers_left.append(x)
            x -= p["a_c"]   
        x_centers_left  = [x - p["g_c"]/2 for x in x_centers_left]

        x_centers_right = []
        x = 0
        for i in range(p["n_right"]):
            x_centers_right.append(x)
            x += p["a_c"]
        x_centers_right = [x + p["g_c"]/2 for x in x_centers_right]

        x_centers = np.concatenate((x_centers_left, x_centers_right))
        y_centers = np.zeros(len(x_centers))

        atom = sm.make_circle(0,0, p["r_c"], layer=1, to_poly=True)
        atoms = sm.GeomGroup()
        for x, y in zip(x_centers, y_centers):
            atoms+=atom.copy().translate(x, y)
        

        cavity = atoms
        return cavity


class Nanobeam(Cavity):
    def initialize(self):
        self.set_name("NANOBEAM")
        self.set_description("A simple nanobeam")

    def parameters(self):
        super().parameters()
        self.addparameter("nanobeam_width", 0.500, "Width of the nanobeam", param_type=float)
        self.addparameter("nanobeam_straight_length", 10.000, "Length of the nanobeam", param_type=float)
        self.addparameter("nanobeam_taper_length", 5.000, "Length of the taper", param_type=float)
        self.addparameter("nanobeam_taper_tip_width", 0.100, "Width of the taper tip", param_type=float)
        self.addparameter("nanobeam_taper_n_points", 100, "Number of points in the taper polygon", param_type=int)  
        self.addparameter("nanobeam_taper_m", 1, "Slope of the taper", param_type=float)
        

    def geom(self):
        p = self.get_params()
        cavity = super().geom()
        cavity.set_layer(2)

        # Central straight part
        rect = sm.make_rect(0,0, p["nanobeam_straight_length"], p["nanobeam_width"], layer=1)
        
        # Taper definition
        x = np.linspace(0, p["nanobeam_taper_length"], p["nanobeam_taper_n_points"])
        w1 = p["nanobeam_width"]
        w2 = p["nanobeam_taper_tip_width"]
        m = p["nanobeam_taper_m"]
        l = p["nanobeam_taper_length"]
        a = (w1-w2)/(l**m)
        w = a*(l-x)**m+w2

        taper = sm.make_tapered_path(x, np.zeros(p["nanobeam_taper_n_points"]), w, 1)
        tapers = taper.copy().translate(p["nanobeam_straight_length"]/2, 0)
        tapers += tapers.copy().mirrorX(0)

        nanobeam  =  rect + tapers + cavity
        
        
        return nanobeam
        
        
        
class NanobeamCoupon(Nanobeam):
    def initialize(self):
        self.set_name("NANOBEAM_COUPON")
        self.set_description("A simple coupon for uTP of a nanobeam")

    def parameters(self):
        super().parameters()
        self.addparameter("theter_width", 0.50, "Width of the taper", param_type=float)
        self.addparameter("theter_length", 5.000, "Length of the taper", param_type=float)
        self.addparameter("theter_both_sides", 1, "Taper on both sides", param_type=int, param_range=[0,1])
        self.addparameter("theter_distance_from_tip", 5.00, "Distance from the tip of the taper", param_type=float)
        self.addparameter("theter_space", 2.00, "Space between the tapers", param_type=float)
        self.addparameter("shrink_um", 0.01, "Shrink in um", param_type=float)

    def geom(self):
        p = self.get_params()
        nanobeam = super().geom()
    
       
        # Theter definition 
        theter_total_length = p["theter_length"] + p["nanobeam_width"]  
        theter = sm.make_rect(0, 0, p["theter_width"], theter_total_length, layer=1, numkey=2)

        theter_space = p["nanobeam_straight_length"] + 2*p["nanobeam_taper_length"] - 2*p["theter_distance_from_tip"]
        theters = sm.GeomGroup()

        x_theters = []
        x = -theter_space/2 
        
        while x < theter_space/2:
            x_theters.append(x)
            x += p["theter_space"]
        x_rest = theter_space/2 - x_theters[-1]
        x_theters = [x + x_rest/2 for x in x_theters]

        cavity = nanobeam.select_layer(2)
        cavity_bb = cavity.bounding_box()
        llx_cavity_bb = cavity_bb.llx
        lrx_cavity_bb = cavity_bb.llx + cavity_bb.width

        x_theters = [x for x in x_theters if not (llx_cavity_bb*0.4 <= x <= lrx_cavity_bb*0.4)]

        for x in x_theters:
            theters += theter.copy().translate(x, 0 )
            
        coupon_width = theter_total_length

        if p["theter_both_sides"]==1:
            theters += theters.copy().mirrorY(0)
            coupon_width *= 2
        else:
            coupon_width += p["nanobeam_width"]+0.3

        theters.set_layer(3)
            
        
        nanobeam_bb = nanobeam.bounding_box()
        llx_nanobeam_bb = nanobeam_bb.llx
        lrx_nanobeam_bb = nanobeam_bb.llx + nanobeam_bb.width

        llx_coupon = llx_nanobeam_bb - 0.1 * nanobeam_bb.width
        lrx_coupon = lrx_nanobeam_bb + 0.1 * nanobeam_bb.width
        coupon_length = lrx_coupon - llx_coupon 

        coupon = sm.make_rect(llx_coupon, theter_total_length, coupon_length, coupon_width, layer=4, numkey=7)

        

        low_current = nanobeam.copy().select_layer(1) + theters.copy().set_layer(1)
        low_current = low_current.poly_outlining(0.300, 1)
        low_current = low_current.boolean_intersection(coupon, 1, 4)

        buffer  = nanobeam.copy().select_layer(1) + theters.copy().set_layer(1)
        buffer = buffer.poly_resize(0.150, 1)
        buffer = buffer.boolean_intersection(coupon, 1, 4)


        high_current = coupon.copy().set_layer(5)
        high_current = high_current.boolean_difference(buffer, 5, 1)

    
        low_current += cavity.copy().set_layer(4)

        low_current.set_layer(1)
        high_current.set_layer(2)

        low_current.poly_resize(-p["shrink_um"], 1)
        high_current.poly_resize(-p["shrink_um"], 2)

       
        id_text = self.make_id_text()  
        annotations = self.write_notes() 
        return low_current + high_current + id_text.translate(0, theter_total_length + 4) + annotations
        
# nanobeam_coupon = NanobeamCoupon.build()
# DeviceInspect(nanobeam_coupon)


import samplemaker.phc as sphc

# Use make_phc_reference to create a reference to the photonic crystal first, call it phc
def make_phc_reference( length, Ny, lattice_constant, radius, mask):
    a = lattice_constant 
    r = radius
    Nx = int(length/a)
    crystal = sphc.Crystal()
    crystal = crystal.triangular_box(Nx = int((Nx-1)//2), Ny = Ny, Nparams=1)
    phc = sphc.make_phc(crystal, scaling=a, cellparams=[r], x0 = 0, y0 = 0) 
    mask.addCell("PHC", phc) 
    return phc

mask = smlay.Mask("dumb")
class PhCNanobeamCoupon(Nanobeam):
    def __init__(self):
        super().__init__()
        self.phc = None

    def initialize(self):
        self.set_name("NANOBEAM_COUPON_PHC")    
        self.set_description("A simple PhC coupon for uTP of a nanobeam")


    def parameters(self):
        super().parameters()
        self.addparameter("theter_width", 0.50, "Width of the taper", param_type=float)
        self.addparameter("theter_length", 1.000, "Length of the taper", param_type=float)
        self.addparameter("theter_separation", 2.000, "Separation between the tapers", param_type=float)
        self.addparameter("pad_edge_x", 40, "Edge of the pad in x", param_type=float)
        self.addparameter("pad_edge_y", 10, "Edge of the pad in y", param_type=float)
        self.addparameter("theter_distance_from_corner", 1.00, "Distance from the corner of the pad", param_type=float) 
        self.addparameter("shrink_um", 0.010, "Shrink factor", param_type=float)
    def set_mirror_phc(self,phc):
        self.phc = phc
        
  
    def geom(self): 
        p = self.get_params()
        nanobeam = super().geom()
        cavity = nanobeam.select_layer(2)
        

        pads = sm.make_rect(0,  0 , p["pad_edge_x"], 2*p["pad_edge_y"]+p["nanobeam_width"], layer=3, numkey=5)
          

        nanobeam_bb = nanobeam.bounding_box()
        llx_nanobeam_bb = nanobeam_bb.llx  
        lrx_nanobeam_bb = nanobeam_bb.llx + nanobeam_bb.width
        llx = llx_nanobeam_bb - 0.5
        lrx = lrx_nanobeam_bb + 0.5
        

        rect = sm.make_rect(0, 0, lrx-llx, nanobeam_bb.height, layer=1, numkey=5)
        

        # theters  in X edge

        theter = sm.make_rect(0, 0, p["theter_width"], p["theter_length"], layer=1, numkey=2)
        x_thetrs = []
        theter_space = p["pad_edge_x"] - 2*p["theter_distance_from_corner"]
        x = -theter_space/2
        while x < theter_space/2:
            x_thetrs.append(x)
            x += p["theter_separation"]
        x_rest = theter_space/2 - x_thetrs[-1]
        x_thetrs = [x + x_rest/2 for x in x_thetrs]

        theters = sm.GeomGroup()
        for x in x_thetrs:
            theters += theter.copy().translate(x, p["pad_edge_y"]+p["nanobeam_width"]/2)
        theters.set_layer(1)
        pads += theters + theters.copy().mirrorY(0) 

        # theters in Y edge
        theter = sm.make_rect(0, 0, p["theter_length"], p["theter_width"], layer=1, numkey=4)

        y_thetrs = []
        theter_space = 2*p["pad_edge_y"] - 2*p["theter_distance_from_corner"] + p["nanobeam_width"]
        y = -theter_space/2
        while y < theter_space/2:
            y_thetrs.append(y)
            y += p["theter_separation"]
        y_rest = theter_space/2 - y_thetrs[-1]
        y_thetrs = [y + y_rest/2 for y in y_thetrs]

        theters = sm.GeomGroup()
        for y in y_thetrs:
            theters += theter.copy().translate(p["pad_edge_x"]/2, y)
        theters.set_layer(1)
        pads += theters + theters.copy().mirrorX(0)
        pads.select_layer(1).set_layer(3)


        # PhC
        
        phc = self.phc if self.phc is not None else make_phc_reference(20, 5, 0.350, 0.080, mask)
        phc_ref = sm.make_sref( 0, 0, "PHC", phc, mag=1, angle=0, mirror= False)
        phc_ref_bb = phc_ref.bounding_box()
        phc_copy =  phc_ref.copy().translate(0,+phc_ref_bb.height/2+p["nanobeam_width"]/2) + phc_ref.copy().translate(0,-phc_ref_bb.height/2 - p["nanobeam_width"]/2)
        phc_copy = phc_copy.flatten()
        phc_copy.set_layer(1)
        phc_copy.all_to_poly(32)



        # nanobeam_low_current & high_current
        nb_low_current = nanobeam.select_layer(1)
        nb_low_current = nb_low_current.poly_outlining(0.300, 1)
        nb_low_current = nb_low_current.boolean_intersection(
            rect.copy().set_layer(1), 1, 1)

        nb_buffer = nanobeam.select_layer(1)
        nb_buffer = nb_buffer.poly_resize(0.150, 1)

        nb_high_current = rect.copy().set_layer(1)
        nb_high_current = nb_high_current.boolean_difference(nb_buffer, 1, 1)


        # pads_low_current & high_current
        pads_bb = pads.bounding_box().toRect().set_layer(1) 

        pads_low_current = pads.copy().set_layer(1).poly_outlining(0.300, 1)
        pads_low_current = pads_low_current.boolean_intersection(pads_bb, 1, 1)
        pads_buffer = pads.copy().set_layer(1).poly_resize(0.150, 1) 
        

        pads_high_current = pads_bb.boolean_difference(pads_buffer, 1, 1)

        high_current = nb_high_current + pads_high_current
        low_current = nb_low_current + pads_low_current + phc_copy + cavity

        high_current.set_layer(2)
        low_current.set_layer(1)

        high_current.poly_resize(-p["shrink_um"],2)
        low_current.poly_resize(-p["shrink_um"],1)


        id_text = self.make_id_text().translate(0, p["pad_edge_y"]+p["nanobeam_width"]/2 + 4)
        annotations = self.write_notes()

        return high_current + low_current + id_text + annotations
    

    
registerDevicesInModule(__name__)

    


       
from itertools import product
group_index = 0
element_row_index = 0
def generate_param_dict(**kwargs):
    # Initialize the dictionary that will store lists for each parameter
    param_dict = {key: [] for key in kwargs}
    
    # Get the cartesian product of all parameter values
    param_combinations = product(*kwargs.values())
    
    # Populate the param_dict with the combinations
    for combination in param_combinations:
        for i, key in enumerate(kwargs):
            param_dict[key].append(combination[i])
    
    return param_dict, len(param_dict[key])



# phc_nb_coupon = PhCNanobeamCoupon()
# DeviceInspect(phc_nb_coupon)

