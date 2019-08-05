# Two component system, with species "A" and "B".
# This follows section 2.3 of:
#
# "A Practical Guide to Stochastic Simulations of
# Reaction-Diffusion Processes" by Erban et al.
#
# results can be compared to Figure 2.3 

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
B = Species('B')
C = Species('C')
D = Species('D')

k1 = 1e-3/unit.sec
k2 = 1e-2/unit.sec
k3 = 1.2/unit.sec
k4 = 1.0/unit.sec

rxns = []
rxns.append(Reaction('AAC',[A],[C],[2],[1],kf=k1))
rxns.append(Reaction('ABD',[A,B],[D],[1,1],[1],kf=k2))
rxns.append(Reaction('birth_A',[],[A],[],[1],kf=k3))
rxns.append(Reaction('birth_B',[],[B],[],[1],kf=k4))

# create a Model
main_compartment = Compartment('main')
main_compartment.add_rxns_to_compartment(rxns)
model = Model(compartments=[main_compartment])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q([0,1,2,3],0)
ode_results = sys.run(100)

fig, ax = plt.subplots(nrows=1,ncols=2)

A_idx = sys.state.index['main']['A']
B_idx = sys.state.index['main']['B']
ax[0].plot(ode_results.t,ode_results.y[A_idx],label='ODE')
ax[1].plot(ode_results.t,ode_results.y[B_idx],label='ODE')

#----
# Now create a GillespieSystem for the same model and run it 10 times
#---

for i in range(5):
    Gillespie_sys = GillespieSystem(flat_model)
    Gillespie_sys.add_reporter(AllReporter(freq=1))
    Gillespie_sys.set_q([0,1,2,3],0)
    tmp = Gillespie_sys.run(100)
    reports = Gillespie_sys.reporters[0].reports()
    t = [r['t'] for r in reports]
    A_idx = Gillespie_sys.state.index['main']['A']
    B_idx = Gillespie_sys.state.index['main']['B']
    A = [r['report'][A_idx] for r in reports]
    B = [r['report'][B_idx] for r in reports]
    ax[0].plot(t,A,label='run {0}'.format(i))
    ax[1].plot(t,B,label='run {0}'.format(i))

# this result should look like Figure 2.1 of Erban et al

ax[0].legend()
ax[1].legend()
plt.show()
