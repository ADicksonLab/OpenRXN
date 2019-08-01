"""States are containers for occupancy values.  States hold a 
number of numpy arrays which are all of the same length.

state.q_val holds a numpy array of quantity values

These are the most important values in the system, they 
show the quantity of material (in moles) of a given species 
in a given compartment.

The rest of the arrays are useful for making selections:

state.species is a numpy (string) array of species IDs
state.compartment is a numpy (string) array of compartment IDs
state.x_pos is a numpy array of x_positions of the compartment centers
state.y_pos is a numpy array of y_positions of the compartment centers
state.z_pos is a numpy array of z_positions of the compartment centers

Note that species values for a given compartment will only be created
if there is either:
1) a Reaction in that compartment which involves that Species
or 2) a connection from that compartment which involves that Species

When initializing the state, an index is constructed to easily 
determine an index given a compartment and a species:

index = state.index[compID][specID]
"""

from openrxn import unit
from openrxn.model import FlatModel
import numpy as np
import pandas as pd

class State(object):
    def __init__(self, model=None, dataframe=None, units=[unit.nanometer]*3):
        """State objects can be initialized using either a 
        FlatModel or a dataframe object.  At minimum, the 
        dataframe needs to have "species" and "compartment" 
        columns."""

        self.index = {}
        self.units = units
        
        if model is not None: 
            assert isinstance(model,FlatModel), "Error! A state object needs a FlatModel to initialize."
            self._init_from_model(model)
        elif dataframe is not None:
            if 'species' not in dataframe.columns or 'compartment' not in dataframe.columns:
                raise ValueError("Error! dataframe must contain columns for 'species' and 'compartment'")
            self._init_from_df(dataframe)        
        
        self.size = len(self.compartment)
        self.q_val = np.zeros((self.size))

    def _init_from_model(self, model):
        
        big_species_list = []
        big_comp_list = []
        big_x_list = []
        big_y_list = []
        big_z_list = []

        running_index = 0        
        for c_tag, c in model.compartments.items():

            self.index[c_tag] = {}
            
            # figure out which species are associated with this compartment
            spec = []
            for other_c, conn in c.connections.items():
                spec += list(conn[1].species_rates.keys())
            for rxn in c.reactions:
                spec += rxn.reactant_IDs
                spec += rxn.product_IDs
            spec_set = set(spec)

            for i,s in enumerate(spec_set):
                self.index[c_tag][s] = running_index + i

            running_index += len(spec_set)

            big_species_list += list(spec_set)
            big_comp_list += [c.ID]*len(spec_set)

            # for x, y and z, average the boundary values
            x = [None,None,None]
            for i in range(len(c.pos)):
                x[i] = 0.5*(c.pos[i][0]+c.pos[i][1]).to(self.units[i]).magnitude
            
            big_x_list += [x[0]]*len(spec_set)
            big_y_list += [x[1]]*len(spec_set)
            big_z_list += [x[2]]*len(spec_set)

        self.species = np.array(big_species_list)
        self.compartment = np.array(big_comp_list)
        self.x_pos = np.array(big_x_list)
        self.y_pos = np.array(big_y_list)
        self.z_pos = np.array(big_z_list)

    def _init_from_df(self, df):

        # assign columns to self arrays
        self.species = np.array(df['species'])
        self.compartment = np.array(df['compartment'])
        if 'x_pos' in df.columns:
            self.x_pos = np.array(df['x_pos'])
        if 'y_pos' in df.columns:
            self.y_pos = np.array(df['y_pos'])
        if 'z_pos' in df.columns:
            self.z_pos = np.array(df['z_pos'])

        # building self.index dictionary
        for i in range(len(df['species'])):
            if df['compartment'][i] not in self.index:
                self.index[df['compartment'][i]] = {}
            self.index[df['compartment'][i]][df['species'][i]] = i

    def to_dataframe(self):
        df = pd.DataFrame()
        
        df['species'] = self.species
        df['compartment'] = self.compartment
        df['q_val'] = self.q_val
        
        if hasattr(self,'x_pos'):
            df['x_pos'] = self.x_pos
        if hasattr(self,'y_pos'):
            df['y_pos'] = self.y_pos
        if hasattr(self,'z_pos'):
            df['z_pos'] = self.z_pos
        
        return df
            
    def to_csv(self, filename):
        df = self.to_dataframe()
        df.to_csv(filename)

    def to_csv_no_q(self, filename):
        df = self.to_dataframe()

        cols = [c for c in df.columns if c != 'q_val']
        df.to_csv(filename,columns=cols)
        
    
