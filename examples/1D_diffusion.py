# One dimensional diffusion system, with one species.
# This follows section 3.2 of:
#
# "A Practical Guide to Stochastic Simulations of
# Reaction-Diffusion Processes" by Erban et al.
#
# Last panel in Figure can be compared to Figure 3.3A
# of Erban et al.
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

# define species and reactions
A = Species('A')

d = 0.16/unit.sec
K = 40   # number of compartments
L = 1*unit.mm
h = L/K

# diffusion constant
D = d*h**2

# create a Model
boundaries = np.linspace(0,L.magnitude,K+1)*L.units
conn = IsotropicConnection({'A' : d},dim=1)
comp_array = CompartmentArray1D('main',boundaries,conn)
model = Model(arrays=[comp_array])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q(np.arange(sys.state.size),0)
sys.set_q([16,17],500)

print("Running ODE...")
ode_results = sys.run(240)

fig, ax = plt.subplots(nrows=2,ncols=2)
ax[0][0].set_title('t=10s')
ax[0][1].set_title('t=60s')
ax[1][0].set_title('t=120s')
ax[1][1].set_title('t=240s')

IDs = []
pos_x = []
for i in range(K):
    c_name = 'main-{0}'.format(i)
    ID = sys.state.index[c_name]['A']
    IDs.append(ID)
    pos_x.append(sys.state.x_pos[ID]/1000000)  # in mm

t1 = np.argwhere(ode_results.t > 10).min()
t2 = np.argwhere(ode_results.t > 60).min()
t3 = np.argwhere(ode_results.t > 120).min()
t4 = ode_results.y.shape[1]-1
ax[0][0].plot(pos_x,ode_results.y[IDs,t1])
ax[0][1].plot(pos_x,ode_results.y[IDs,t2])
ax[1][0].plot(pos_x,ode_results.y[IDs,t3])
ax[1][1].plot(pos_x,ode_results.y[IDs,t4],label='ODE')

#----
# Now create a GillespieSystem for the same model
#---

for i in range(5):
    print(f"Running Gillespie system (run {i+1} of 5)")
    Gillespie_sys = GillespieSystem(flat_model)
    Gillespie_sys.add_reporter(AllReporter(freq=1))
    Gillespie_sys.set_q(np.arange(sys.state.size),0)
    Gillespie_sys.set_q([16,17],500)
    tmp = Gillespie_sys.run(240)
    reports = Gillespie_sys.reporters[0].reports()
    ax[0][0].plot(pos_x,reports[10]['report'])
    ax[0][1].plot(pos_x,reports[60]['report'])
    ax[1][0].plot(pos_x,reports[120]['report'])
    ax[1][1].plot(pos_x,reports[239]['report'],label='run {0}'.format(i))

ax[1][1].legend()
plt.tight_layout()
outname = "1D_diffusion.pdf"
plt.savefig(outname)

print("Output saved as ",outname)
