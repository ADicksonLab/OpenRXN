# One dimensional diffusion system, with four species, A
# B, C and D, which interact as in Gillespie_AB_system.py
#
#            k1             k2
#      A + A -> C     A + B -> D
#
# A is created only in 0 < x < 9L/10, and
# B is created only in 2L/5 < x < L
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

# define species and reactions
A = Species('A')
B = Species('B')
C = Species('C')
D = Species('D')

# h is the "volume" of the compartment here
# second order rate coefficients should have units of
# M^-1 s^-1, or L^d / s, where L is units of length and
# d is the dimensionality of the system (here, 1)

k1 = 1e-3/unit.sec*h
k2 = 1e-2/unit.sec*h
k3 = 1.2/unit.sec
k4 = 1.0/unit.sec

rxns = []
rxns.append(Reaction('AAC',[A],[C],[2],[1],kf=k1))
rxns.append(Reaction('ABD',[A,B],[D],[1,1],[1],kf=k2))
rxns.append(Reaction('birth_A',[],[A],[],[1],kf=k3))
rxns.append(Reaction('birth_B',[],[B],[],[1],kf=k4))

# create a Model
boundaries = np.linspace(0,L.magnitude,K+1)*L.units
conn = IsotropicConnection({'A' : d*h, 'B' : d*h},dim=1)
comp_array = CompartmentArray1D('main',boundaries,conn)
comp_array.add_rxns_to_array([rxns[0],rxns[1]])
for c in comp_array.compartments.values():
    if c.pos[0][1] <= 9*L/10:
        c.add_rxn_to_compartment(rxns[2])
    if c.pos[0][1] > 2*L/5:
        c.add_rxn_to_compartment(rxns[3])

model = Model(arrays=[comp_array])
flat_model = model.flatten()

# create a system
sys = ODESystem(flat_model)

# set initial concentrations
sys.set_q(np.arange(sys.state.size),0)

ode_results = sys.run(1800)

fig, ax = plt.subplots(nrows=1,ncols=2)
ax[0].set_ylabel('Number of A molecules')
ax[1].set_ylabel('Number of B molecules')
ax[0].set_xlabel('x (mm)')
ax[1].set_xlabel('x (mm)')
plt.show()

IDs_A = []
IDs_B = []
pos_x = []
for i in range(K):
    c_name = 'main-{0}'.format(i)
    IDs_A.append(sys.state.index[c_name]['A'])
    IDs_B.append(sys.state.index[c_name]['B'])
pos_x = sys.state.x_pos[IDs_A]/1000000  # in mm

t4 = ode_results.y.shape[1]-1
ax[0].plot(pos_x,ode_results.y[IDs_A,t4],label='ODE')
ax[1].plot(pos_x,ode_results.y[IDs_B,t4],label='ODE')

#----
# Now create a GillespieSystem for the same model
#---

for i in range(1):
    Gillespie_sys = GillespieSystem(flat_model)
    Gillespie_sys.add_reporter(AllReporter(freq=100))
    Gillespie_sys.set_q(np.arange(sys.state.size),0)
    tmp = Gillespie_sys.run(1800)
    IDs_A = []
    IDs_B = []
    for j in range(K):
        c_name = 'main-{0}'.format(j)
        IDs_A.append(sys.state.index[c_name]['A'])
        IDs_B.append(sys.state.index[c_name]['B'])
    reports = Gillespie_sys.reporters[0].reports()
    ax[0].plot(pos_x,reports[-1]['report'][IDs_A],label='run {0}'.format(i))
    ax[1].plot(pos_x,reports[-1]['report'][IDs_B],label='run {0}'.format(i))

ax[0].legend()
ax[1].legend()
plt.show()
