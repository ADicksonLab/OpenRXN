# Very, very simple example showing the degradation of
# species "A" in a single compartment.  This follows
# section 2.1 of:
#
# "A Practical Guide to Stochastic Simulations of
# Reaction-Diffusion Processes" by Erban et al.
#
# results can be compared to Figure 2.1

from openrxn.systems.ODESystem import ODESystem
from openrxn.reactions import Reaction, Species
from openrxn.model import Model
from openrxn.compartments.compartment import Compartment
from openrxn import unit

# define species and reactions
A = Species('A')
k = 0.1/unit.sec
degradation = Reaction('degradation',[A],[],[1],[],kf=k)

# create a Model
main_compartment = Compartment('main')
main_compartment.add_rxn_to_compartment(degradation)
model = Model(compartments=[main_compartment])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q([0],20)
results = sys.run(30)

# this result should look like Figure 2.1 of Erban et al
import matplotlib.pyplot as plt
plt.plot(results.t,results.y[0])

plt.show()
