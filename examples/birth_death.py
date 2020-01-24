# Simple example showing the birth and death of
# species "A" in a single compartment.  This follows
# section 2.2 of:
#
# "A Practical Guide to Stochastic Simulations of
# Reaction-Diffusion Processes" by Erban et al.
#
# results can be compared to Figure 2.2

from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.GillespieSystem import GillespieSystem
from openrxn.reactions import Reaction, Species
from openrxn.model import Model
from openrxn.compartments.compartment import Compartment
from openrxn import unit
from openrxn.systems.reporters import AllReporter, SumReporter, SelectionReporter
import matplotlib.pyplot as plt

# define species and reactions
A = Species('A')
k1 = 0.1/unit.sec
k2 = 1.0/unit.sec
birth_and_death = Reaction('birth_and_death',[A],[],[1],[],kf=k1,kr=k2)

# create a Model
main_compartment = Compartment('main')
main_compartment.add_rxn_to_compartment(birth_and_death)
model = Model(compartments=[main_compartment])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q([0],0)
ode_results = sys.run(100)
plt.plot(ode_results.t,ode_results.y[0],label='ODE')

#----
# Now create a GillespieSystem for the same model and run it 10 times
#---

for i in range(10):
    Gillespie_sys = GillespieSystem(flat_model)
    Gillespie_sys.add_reporter(AllReporter(freq=1))
    Gillespie_sys.set_q([0],0)
    tmp = Gillespie_sys.run(100)
    reports = Gillespie_sys.reporters[0].reports()
    t = [r['t'] for r in reports]
    y = [r['report'] for r in reports]
    plt.plot(t,y,label='run {0}'.format(i))

# this result should look like Figure 2.1 of Erban et al
plt.show()



