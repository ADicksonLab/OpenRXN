# One dimensional diffusion system, with one species, A
# which is degraded everywhere, but created only in
# compartments with x <= L/5.
#
# This follows section 4.1 of:
#
# "A Practical Guide to Stochastic Simulations of
# Reaction-Diffusion Processes" by Erban et al.
#
# The last two panels in the Figure can be compared to
# Figure 4.1 of Erban et al.
#

from openrxn.systems.ODESystem import ODESystem
from openrxn.systems.GillespieSystem import GillespieSystem
from openrxn.reactions import Reaction, Species
from openrxn.model import Model
from openrxn.compartments.arrays import CompartmentArray1D
from openrxn import unit
from openrxn.systems.reporters import AllReporter, SumReporter, SelectionReporter
from openrxn.connections import IsotropicConnection

import matplotlib.pyplot as plt
import numpy as np

d = 0.16/unit.sec
K = 40   # number of compartments
L = 1*unit.mm
h = L/K

# diffusion constant
D = d*h**2

# define species and reactions
A = Species('A')
kp = 0.012/(unit.micrometer*unit.sec)
rxns = []
rxns.append(Reaction('deg',[A],[],[1],[],kf=1e-3/unit.sec))
rxns.append(Reaction('syn',[],[A],[],[1],kf=kp*h))

# create a Model
boundaries = np.linspace(0,L.magnitude,K+1)*L.units
conn = IsotropicConnection({'A' : d},dim=1)
comp_array = CompartmentArray1D('main',boundaries,conn)
comp_array.add_rxn_to_array(rxns[0])
for c in comp_array.compartments.values():
    if c.pos[0][1] <= L/5:
        c.add_rxn_to_compartment(rxns[1])

model = Model(arrays=[comp_array])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q(np.arange(sys.state.size),0)

ode_results = sys.run(1800)

fig, ax = plt.subplots(nrows=2,ncols=2)
ax[0][0].set_title('t=1 min')
ax[0][1].set_title('t=5 min')
ax[1][0].set_title('t=10 min')
ax[1][1].set_title('t=30 min')

IDs = []
pos_x = []
for i in range(K):
    c_name = 'main-{0}'.format(i)
    ID = sys.state.index[c_name]['A']
    IDs.append(ID)
    pos_x.append(sys.state.x_pos[ID]/1000000)  # in mm

t1 = np.argwhere(ode_results.t > 60).min()
t2 = np.argwhere(ode_results.t > 300).min()
t3 = np.argwhere(ode_results.t > 600).min()
t4 = ode_results.y.shape[1]-1
ax[0][0].plot(pos_x,ode_results.y[IDs,t1])
ax[0][1].plot(pos_x,ode_results.y[IDs,t2])
ax[1][0].plot(pos_x,ode_results.y[IDs,t3])
ax[1][1].plot(pos_x,ode_results.y[IDs,t4],label='ODE')

#----
# Now create a GillespieSystem for the same model
#---

for i in range(1):
    Gillespie_sys = GillespieSystem(flat_model)
    Gillespie_sys.add_reporter(AllReporter(freq=10))
    Gillespie_sys.set_q(np.arange(sys.state.size),0)
    tmp = Gillespie_sys.run(1800)
    reports = Gillespie_sys.reporters[0].reports()
    ax[0][0].plot(pos_x,reports[5]['report'])
    ax[0][1].plot(pos_x,reports[29]['report'])
    ax[1][0].plot(pos_x,reports[59]['report'])
    ax[1][1].plot(pos_x,reports[179]['report'],label='run {0}'.format(i))

ax[1][1].legend()
plt.show()
