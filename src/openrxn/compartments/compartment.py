"""Compartments hold a set of reactions and govern the transport
of material between compartments through connections."""

from openrxn.compartments.ID import make_ID
from openrxn import unit
import logging

class Compartment(object):
    """Compartments are initialized with an ID, which can be a string, an int
    or a tuple.

    pos is an optional argument for the position of the compartment,
    which can be useful to determine flow rates and to construct spatial maps.
    pos is array-like, where the length of the array denotes the 
    dimensionality of the space.

    self.reactions is a list of reactions that occur in this compartment

    self.connections is a list of compartment IDs that this compartment
    is connected to.

    Upon initialization, both of these lists are empty.
    """
    
    def __init__(self, ID, pos=None, array_ID=None, volume=None):
        self.ID = ID
        self._rxn_ids = set()
        self.reactions = []
        self.connections = {}
        self.pos = [] if pos is None else list(pos)
        self.array_ID = array_ID
        self.volume = volume
        
    def add_rxn_to_compartment(self, rxn):
        """Adds a reaction to a compartment."""
        if rxn.ID in self._rxn_ids:
            logging.warning("Reaction {0} already in compartment {1}".format(rxn.ID,self.ID))
        else:
            self.reactions.append(rxn)
            self._rxn_ids.add(rxn.ID)

    def add_rxns_to_compartment(self, rxns):
        """Adds a list of reactions to a compartment."""
        for rxn in rxns:
            self.add_rxn_to_compartment(rxn)
            
    def show_all_rxns(self):
        """Returns a list of rxn strings with all reactions in the compartment"""
        rxn_strings = []
        for r in self.reactions:
            rxn_strings.append(r.display())
            
        return rxn_strings
        
    def connect(self, other_compartment, conn_type, warn_overwrite=True):
        """Make a connection from this compartment to another one
        using the conn_type connection type."""

        conn_tag = make_ID(other_compartment.array_ID,other_compartment.ID)
        if conn_tag in self.connections and warn_overwrite:
            self_tag = make_ID(self.array_ID,self.ID)
            logging.warning("Overwriting connection between {0} and {1}".format(self_tag,conn_tag))

        self.connections[conn_tag] = (other_compartment, conn_type)

    def remove_connection(self, other_compartment):
        """Remove the connection with the other_compartment"""

        conn_tag = make_ID(other_compartment.array_ID,other_compartment.ID)
        if conn_tag not in self.connections:
            self_tag = make_ID(self.array_ID,self.ID)
            logging.warning("Connection to remove between {0} and {1} does not exist".format(self_tag,conn_tag))

        val = self.connections.pop(conn_tag)

    def copy(self,ID=None,delete_array_ID=False):
        """Returns an identical copy of this compartment."""

        if ID == None:
            newID = self.ID
        else:
            newID = ID
        if delete_array_ID:
            new_aID = None
        else:
            new_aID = self.array_ID

        if hasattr(self,'surface_area'):
            new_comp = type(self)(newID, pos=self.pos, array_ID=new_aID, surface_area=self.surface_area)
        else:
            new_comp = type(self)(newID, pos=self.pos, array_ID=new_aID)
            
        new_comp.volume = self.volume
        new_comp.connections = dict(self.connections)
        new_comp.reactions = list(self.reactions)
        new_comp._rxn_ids = {r.ID for r in new_comp.reactions}

        return new_comp

    def __repr__(self):
        rxn_strings = self.show_all_rxns()
        vol_str = str(self.volume) if self.volume is not None else "None"
        pos_str = tuple(self.pos) if self.pos else None

        return (
            f"Compartment(ID={self.ID!r}, "
            f"volume={vol_str}, pos={pos_str}, "
            f"n_rxns={len(self.reactions)}, n_conns={len(self.connections)})"
        )

class Compartment1D(Compartment):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        assert len(self.pos)==1, "Error! position must be of length 1 for Compartment1D objects"
        self.volume = (self.pos[0][1]-self.pos[0][0])
        self.volume.to(unit.nm)

class Compartment2D(Compartment):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        assert len(self.pos)==2, "Error! position must be of length 2 for Compartment2D objects"
        self.volume = (self.pos[0][1]-self.pos[0][0])*(self.pos[1][1]-self.pos[1][0])
        self.volume.to(unit.nm**2)

class Compartment3D(Compartment):
    """
    surface_area : dict
    An optional argument to describe the surface area of the faces
    'xy', 'yz', 'xz'.  Used by CompartmentArray3D for FicksConnections.
    """

    def __init__(self, *args, surface_area=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.surface_area = surface_area
        assert len(self.pos)==3, "Error! position must be of length 3 for Compartment3D objects"
        self.volume = (self.pos[0][1]-self.pos[0][0])*(self.pos[1][1]-self.pos[1][0])*(self.pos[2][1]-self.pos[2][0])
        self.volume.to(unit.nm**3)

class Reservoir(Compartment):
    """
    A compartment with a fixed concentration of different species.
    """

    def __init__(self, *args, concs={}, conc_funcs={}, **kwargs):
        super().__init__(*args, **kwargs)

        self.conc_funcs = {}

        # add keys from concs into self.conc_funcs
        for key in concs.keys():
            # convert concentrations to number of molecules per
            # nanometer**3
            concs[key].ito(1/unit.nm**3)

            ## define a concentration function and attach it
            def tmp_conc_func(t):
                return concs[key].ito(1/unit.nm**3)
            
            self.conc_funcs[key] = tmp_conc_func

        # add keys from conc_funcs into self.conc_funcs
        for key in conc_funcs.keys():
            if key in concs:
                raise ValueError("Error! Same quantity passed in both concs and conc_funcs!")
            if not callable(conc_funcs[key]):
                raise ValueError("Error! conc_funcs[{0}] is not callable!".format(key))
            self.conc_funcs[key] = conc_funcs[key]

    def copy(self,ID=None,delete_array_ID=False):
        """Returns an identical copy of this compartment."""

        if ID == None:
            newID = self.ID
        else:
            newID = ID
        if delete_array_ID:
            new_aID = None
        else:
            new_aID = self.array_ID

        new_comp = Reservoir(newID, array_ID=new_aID, conc_funcs=self.conc_funcs)
            
        return new_comp
