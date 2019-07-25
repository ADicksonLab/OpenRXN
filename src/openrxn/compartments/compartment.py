"""Compartments hold a set of reactions and govern the transport
of material between compartments through connections."""

from openrxn.compartments.ID import makeID

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

    surface_area : dict
    An optional argument to describe the surface area of the faces
    'xy', 'yz', 'xz'.  Used by CompartmentArray3D for FicksConnections.
    """
    
    def __init__(self, ID, pos=None, array_ID=None, surface_area=None):
        self.ID = ID
        self.reactions = []
        self.connections = {}
        self.pos = pos
        self.array_ID = array_ID
        self.surface_area = surface_area
        
    def add_rxn_to_compartment(self, rxn):
        """Adds a reaction to a compartment."""
        if rxn.ID in [r.ID for r in self.reactions]:
            logging.warn("Reaction {0} already in compartment {1}".format(rxn.ID,self.ID))
        else:
            self.reactions.append(rxn)

    def show_all_rxns(self):
        """Returns a list of rxn strings with all reactions in the compartment"""
        rxn_strings = []
        for r in self.reactions:
            rxn_strings.append(r.display())
            
        return rxn_strings
        
    def connect(self, other_compartment, conn_type, warn_overwrite=True):
        """Make a connection from this compartment to another one
        using the conn_type connection type."""

        conn_tag = makeID(other_compartment.array_ID,other_compartment.ID)
        if conn_tag in self.connections and warn_overwrite:
            self_tag = makeID(self.array_ID,self.ID)
            logging.warn("Warning: overwriting connection between {0} and {1}".format(self_tag,conn_tag))

        self.connections[conn_tag] = (other_compartment, conn_type)

    def remove_connection(self, other_compartment):
        """Remove the connection with the other_compartment"""

        conn_tag = makeID(other_compartment.array_ID,other_compartment.ID)
        if conn_tag not in self.connections:
            self_tag = makeID(self.array_ID,self.ID)
            logging.warn("Warning: connection to remove between {0} and {1} does not exist".format(self_tag,conn_tag))

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
        
        new_comp = Compartment(newID, pos=self.pos, array_ID=new_aID, surface_area=self.surface_area)
        new_comp.connections = self.connections
        new_comp.reactions = self.reactions

        return new_comp

            
