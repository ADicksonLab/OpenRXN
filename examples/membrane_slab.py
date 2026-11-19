# This example sets up a membrane slab with some receptors,
# with solvent above, and a top reservoir compartment

from openrxn.systems.ODESystem import ODESystem
from openrxn.reactions import Reaction, Species
from openrxn.compartments.arrays import CompartmentArray3D
from openrxn.connections import IsotropicConnection, AnisotropicConnection, FicksConnection
from openrxn.model import Model
from openrxn import unit
from openrxn.systems.reporters import AllReporter, SumReporter, SelectionReporter

import os
import numpy as np
import networkx as nx

# define species and reactions

drug = Species('drug')
receptor = Species('receptor')
dr_complex  = Species('complex')

kon = 1e6/(unit.mol*unit.sec/unit.liter)
koff = 0.1/unit.sec

binding = Reaction('binding',[drug,receptor],[dr_complex],[1,1],[1],kf=kon,kr=koff)

# define connections and compartments

in_slab = FicksConnection({'drug' : 1e-8*unit.cm**2/unit.sec})

x_pos = np.linspace(-50,50,11)*unit.nanometer
y_pos = np.linspace(-50,50,11)*unit.nanometer
z_pos1 = np.array([-1,0])*unit.nanometer
lower_slab = CompartmentArray3D('lower_slab',x_pos,y_pos,z_pos1,in_slab,periodic=[True,True,False])

z_pos2 = np.array([0,1])*unit.nanometer
upper_slab = CompartmentArray3D('upper_slab',x_pos,y_pos,z_pos2,in_slab,periodic=[True,True,False])
upper_slab.add_rxn_to_array(binding)

between_slab = IsotropicConnection({'drug' : 1e-5/unit.sec})
lower_slab.join3D(upper_slab,between_slab,append_side='z+')

# define compartment for bulk
bulk_to_bulk = FicksConnection({'drug' : 1e-5*unit.cm**2/unit.sec})
z_pos3 = np.linspace(1,25,23)*unit.nanometer
bulk = CompartmentArray3D('bulk',x_pos,y_pos,z_pos3,bulk_to_bulk,periodic=[True,True,False])

# (k_out, k_in) is read from bulk's side: bulk -> slab at 1e-5 /s, slab -> bulk at 1e-1 /s;
# join3D gives the slab side the reversed rates
slab_to_bulk = AnisotropicConnection({'drug' : (1e-5/unit.sec, 1e-1/unit.sec)})
bulk.join3D(upper_slab,slab_to_bulk,append_side='z-')

# create model
my_model = Model([lower_slab,upper_slab,bulk])
my_flat_model = my_model.flatten()

# export a graph
graph = my_flat_model.to_graph()
nx.readwrite.gexf.write_gexf(graph,os.path.join(os.path.dirname(os.path.abspath(__file__)),'membrane.gexf'))

# make a system
system = ODESystem(my_flat_model)

# set initial concentrations (state positions are compartment centers, as floats in nm)
select_drug_top_layer = np.where(np.logical_and(
    system.state.z_pos > 23, system.state.species == drug.ID))[0]
system.set_q(select_drug_top_layer,1e-2*unit.mol/unit.L)

select_receptor_top_slab = np.where(np.logical_and(
    np.isclose(system.state.z_pos, 0.5), system.state.species == receptor.ID))[0]
system.set_q(select_receptor_top_slab,1e-2*unit.mol/unit.L)

# integrate
select_drug_all = np.where(system.state.species == drug.ID)[0]
drug_reporter = SelectionReporter(select_drug_all,freq=1e-5)

select_complex_all = np.where(system.state.species == dr_complex.ID)[0]
complex_all_reporter = SumReporter(select_complex_all,freq=1e-5)

system.add_reporters([drug_reporter,complex_all_reporter])
results = system.run(1e-4)

# drug is conserved: free drug + complex stays at its initial amount
n_drug0 = drug_reporter.reports()[0]['report'].sum()
n_drug_end = results.y[select_drug_all, -1].sum() + results.y[select_complex_all, -1].sum()
print("compartments: {0}, states: {1}".format(len(my_flat_model.compartments), system.state.size))
print("drug molecules (free + bound): t=0 {0:.6g}, t=1e-4 s {1:.6g}".format(n_drug0, n_drug_end))
print("complex at report times:", ["{0:.3g}".format(r['report']) for r in complex_all_reporter.reports()])

