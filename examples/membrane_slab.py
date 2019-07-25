# This example sets up a membrane slab with some receptors,
# with solvent above, and a top reservoir compartment

from openrxn.reactions import Reaction, Species
from openrxn.compartments.arrays import CompartmentArray3D
from openrxn.connections import IsotropicConnection, AnisotropicConnection, FicksConnection
from openrxn.model import Model

import pint
unit = pint.UnitRegistry()

import numpy as np

# define species and reactions

drug = Species('drug')
receptor = Species('receptor')
dr_complex  = Species('complex')

kon = 1e6/(unit.molar*unit.sec)
koff = 0.1/unit.sec

binding = Reaction('binding',[drug,receptor],[dr_complex],[1,1],[1],kf=kon,kr=koff)

# define connections and compartments

in_slab = FicksConnection({'drug' : 1e-8*unit.cm**2/unit.sec})

x_pos = np.linspace(-50,50,101)*unit.nanometer
y_pos = np.linspace(-50,50,101)*unit.nanometer
z_pos1 = np.array([-1,0])*unit.nanometer
lower_slab = CompartmentArray3D('lower_slab',x_pos,y_pos,z_pos1,in_slab,periodic=[True,True,False])

z_pos2 = np.array([0,1])*unit.nanometer
upper_slab = CompartmentArray3D('upper_slab',x_pos,y_pos,z_pos2,in_slab,periodic=[True,True,False])

between_slab = IsotropicConnection({'drug' : 1e-5/unit.sec})
lower_slab.join3D(upper_slab,between_slab,append_side='z+')

# define compartment for bulk
bulk_to_bulk = FicksConnection({'drug' : 1e-5*unit.cm**2/unit.sec})
z_pos3 = np.linspace(1,25,23)*unit.nanometer
bulk = CompartmentArray3D('bulk',x_pos,y_pos,z_pos3,bulk_to_bulk,periodic=[True,True,False])

slab_to_bulk = AnisotropicConnection({'drug' : (1e-5/unit.sec, 1e-1/unit.sec)})
bulk.join3D(upper_slab,slab_to_bulk,append_side='z-')

# create model
my_model = Model([lower_slab,upper_slab,bulk])
my_flat_model = my_model.flatten()


