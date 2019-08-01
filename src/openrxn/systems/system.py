"""Systems are fully "realized" models, with a concentration
values specified for each compartment at a given time point
in system.state .
A system is moved forward in time with the system.propagate()
function.

In detail, all systems use a list of functions (system.dqdt), 
where system.dqdt[i] returns the rate of change of quantity i.

The system need not store the values of every species in every 
compartment.  Reporter functions can be attached to system objects
which can return the sum or average the concentrations in different
sets of compartments.  Reporters can also be configured to return 
ALL of the data and this is enabled by default.

After running, results (a.k.a. the reports from the reporters) are 
stored in system.results."""

from openrxn import unit
from openrxn.systems.state import State

import numpy as np
import logging

EPSILON = 1e-8

class System(object):

    def __init__(self, flatmodel, init_state=None, reporters=[]):
        """Systems must be initialized with FlatModel objects.
        initial states can be specified in the init_state argument,
        but care must be taken to ensure that this is compatible
        with the Model.

        It is strongly recommended to let the System initialize its
        own State, and then set its initial values using a set of
        selections and assignments.

        e.g. 
        s = System(my_flatmodel)
        species_a_bottom_layer = np.where(np.logical_and(
                                    s.state.z_pos < 1, s.state.species == a.ID))
        s.state.set_q(species_a_bottom_layer, 1 * ureg.mol)
        """

        self.model = flatmodel
        
        if init_state != None:
            self.state = init_state
        else:
            self.state = State(model=self.model)

        self.reporters = []
        self.reporters += reporters

    def add_reporter(self,reporter):
        self.reporters.append(reporter)

    def add_reporters(self,reporters):
        self.reporters += reporters

    def run(self,total_time,**kwargs):
        """
        Runs the system forward in time using the system-specific
        self.propagate function,

        Takes the total time of integration as an argument.  Other
        keyword arguments (kwargs) are passed to the self.propagate
        function.

        Note:  this function assumes that the initial state vector (q_val)
        has already been set, and that reporters have already been 
        defined and attached to the system.
        """

        report_freqs = [r.freq for r in self.reporters]

        # add endpoints as default, even if reporters aren't present
        checkpoints = [0,total_time]
        for freq in report_freqs:
            n = int(total_time/freq) + 1
            checkpoints += [freq*i for i in range(n)]

        checkpoints = list(set(checkpoints))
        checkpoints.sort()
        
        for i in range(len(checkpoints)-1):
            init_t = checkpoints[i]
            final_t = checkpoints[i+1]
            
            result = self.propagate((init_t,final_t),**kwargs)
            if 'final_t' in result:
                checkpoints[i+1] = result['final_t']

            logging.info("Reached checkpoint: t = {0}".format(checkpoints[i+1]))
            
            for r in self.reporters:
                # check whether final_t is a multiple of r.freq
                if final_t/r.freq - int(final_t/r.freq) < EPSILON:
                    r.report(checkpoints[i+1], self.state.q_val)

        return result
        
    def propagate(self,**kwargs):
        raise NotImplementedError
        
        
