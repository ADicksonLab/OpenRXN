"""Compartment arrays are groups of compartments that make
it easier to connect and manipulate large groups."""

from openrxn.connections import Connection
from openrxn.compartments.compartment import Compartment1D, Compartment2D, Compartment3D

class CompartmentArray(object):
    """Base class for compartment arrays."""

    def to_graph(self):
        """Returns a graph, with compartments as nodes and
        connections as edges."""

    def add_rxn_to_array(self, rxn):
        """Adds a reaction to each compartment in the array."""
        for c in self.compartments.values():
            c.add_rxn_to_compartment(rxn)

    def add_rxns_to_array(self, rxns):
        """Adds each reaction in rxns to each compartment in the array."""
        for r in rxns:
            self.add_rxn_to_array(r)

    def change_all_intra_connection_type(self, new_ctype):
        """Change the connection type between the compartments,
        overwriting them with new_ctype, which must be a Connection 
        type. 

        This only affects connections within the array,
        not connections between arrays."""

        for c in self.compartments.values():
            for n,c_type in c.connections.values():

                # check if both compartments are in same array
                if n.array_ID == self.array_ID:
                    c.connect(n,new_ctype,warn_overwrite=False)    

    def change_all_inter_connection_type(self, other_array, new_ctype):
        """Change the connection type between the compartments in this
        array and those in "other_array", overwriting them with 
        new_ctype, which must be a Connection type. 

        This only affects connections between the arrays,
        not connections within arrays."""

        for c in self.compartments.values():
            for n,c_type in c.connections.values():

                # check if other compartment is in other_array
                if n.array_ID == other_array.array_ID:
                    c.connect(n,new_ctype,warn_overwrite=False)    

                    
class CompartmentArray1D(CompartmentArray):
    """Uses a 1D array of compartments, which are connected
    to each other in sequence.  Their positions must be 
    specified upon initialization and are added as args
    to the compartments.

    the positions list includes the endpoints, so the number 
    of compartments is equal to len(positions) - 1

    conn_type is the Connection to be used to connect the 
    compartments in the array."""

    def __init__(self, array_ID, positions, conn_type, periodic=False):

        self.array_ID = array_ID
        self.box_len = [positions[-1]-positions[0]]
        self.periodic = periodic
        
        assert isinstance(conn_type,Connection), "conn_type must be of type Connection"

        # set number of compartments
        self.n_compartments = len(positions)-1

        # initialize compartment dictionary
        self.compartments = {}
        for i in range(len(positions)-1):
            self.compartments[(i)] = Compartment1D((i), pos=[(positions[i],positions[i+1])], array_ID=self.array_ID)

        # add connections
        for i in range(self.n_compartments):
            if i > 0:
                self.compartments[(i)].connect(self.compartments[(i-1)],conn_type)
            if i < self.n_compartments-1:
                self.compartments[(i)].connect(self.compartments[(i+1)],conn_type)
        if periodic:
            self.compartments[(0)].connect(self.compartments[(self.n_compartments-1)],conn_type)
            self.compartments[(self.n_compartments-1)].connect(self.compartments[(0)],conn_type)

    def stack(self,other_array,conn_type):
        """Stacks another 1D compartment array on top of this
        one, making connections between the two groups

        conn_type is the Connection to be used to connect the 
        compartments between the two arrays."""

        assert isinstance(other_array, CompartmentArray1D), "type of other_array must be CompartmentArray1D"

        assert self.n_compartments == other_array.n_compartments, "other_array must be the same size as this one"

        assert isinstance(conn_type, Connection), "conn_type must be of type Connection"

        # add connections between compartments
        for i in range(self.n_compartments):
            self.compartments[(i)].connect(other_array.compartments[(i)],conn_type)
            other_array.compartments[(i)].connect(self.compartments[(i)],conn_type)

class CompartmentArray2D(CompartmentArray):
    """Uses a 2D array of compartments, which are connected
    to each other in a grid.  Their positions must be 
    specified upon initialization and are added as args
    to the compartments.

    x_pos is a list of x_positions in the grid
    y_pos is a list of y_positions in the grid

    both x_pos and y_pos include endpoints, so
    the total number of compartments in the grid is equal 
    to (len(x_pos)-1) * (len(y_pos)-1)

    conn_type is the Connection to be used to connect the 
    compartments in the array.

    periodic is a list with two boolean elements controlling 
    periodicity in the x and y dimensions respectively.
    """

    def __init__(self, array_ID, x_pos, y_pos, conn_type, periodic=[False,False]):

        self.array_ID = array_ID
        
        assert isinstance(conn_type, Connection), "conn_type must be of type Connection"

        # set number of compartments
        self.nx = len(x_pos)-1
        self.ny = len(y_pos)-1
        self.n_compartments = self.nx*self.ny
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.box_len = [x_pos[-1]-x_pos[0],y_pos[-1]-y_pos[0]]
        self.periodic = periodic

        # initialize compartment dictionary
        self.compartments = {}
        for i in range(self.nx-1):
            posx = (x_pos[i],x_pos[i+1])
            for j in range(self.ny-1):
                posy = (y_pos[j],y_pos[j+1])
                self.compartments[(i,j)] = Compartment2D((i,j), pos=[posx,posy], array_ID=self.array_ID)

        # add connections
        for i in range(self.nx):
            for j in range(self.ny):
                if i > 0:
                    # make connection to the left
                    self.compartments[(i,j)].connect(self.compartments[(i-1,j)],conn_type)
                if i < self.nx-1:
                    # make connection to the right
                    self.compartments[(i,j)].connect(self.compartments[(i+1,j)],conn_type)
                if j > 0:
                    # make connection to the bottom
                    self.compartments[(i,j)].connect(self.compartments[(i,j-1)],conn_type)
                if j < self.ny-1:
                    # make connection to the top
                    self.compartments[(i,j)].connect(self.compartments[(i,j+1)],conn_type)

        if periodic[0]:
            for j in range(self.ny):
                self.compartments[(0,j)].connect(self.compartments[(self.nx-1,j)],conn_type)
                self.compartments[(self.nx-1,j)].connect(self.compartments[(0,j)],conn_type)

        if periodic[1]:
            for i in range(self.nx):
                self.compartments[(i,0)].connect(self.compartments[(i,self.ny-1)],conn_type)
                self.compartments[(i,self.ny-1)].connect(self.compartments[(i,0)],conn_type)
        
    def stack(self,other_array,conn_type):
        """Stacks another 2D compartment array on top of this
        one, making connections between the two groups

        conn_type is the Connection to be used to connect the 
        compartments between the two arrays."""

        assert isinstance(other_array, CompartmentArray2D), "type of other_array must be CompartmentArray2D"

        assert len(self.x_pos) == len(other_array.x_pos), "other_array must be the same size as this one (x)"

        assert len(self.y_pos) == len(other_array.y_pos), "other_array must be the same size as this one (y)"
        
        assert isinstance(conn_type, Connection), "conn_type must be of type Connection"

        # add connections between compartments
        for i in range(self.nx):
            for j in range(self.ny):
                self.compartments[(i,j)].connect(other_array.compartments[(i,j)],conn_type)
                other_array.compartments[(i,j)].connect(self.compartments[(i,j)],conn_type)
        
class CompartmentArray3D(CompartmentArray):
    """Uses a 3D array of cubic compartments, which are connected
    to each other in a grid.  Their positions must be 
    specified upon initialization and are added as args
    to the compartments.

    x_pos is a list of x_positions in the grid 
    y_pos is a list of y_positions in the grid
    z_pos is a list of z_positions in the grid

    all of the pos arrays include endpoints
    (e.g. the x coordinates of the ith cell goes from 
    x_pos[i] to x_pos[i+1])

    the total number of compartments in the grid is equal 
    to len(x_pos)-1 * len(y_pos)-1 * len(z_pos)-1

    conn_type is the Connection to be used to connect the 
    compartments in the array.

    periodic is a list with three boolean elements controlling 
    periodicity in the x, y and z dimensions respectively.
    """

    def __init__(self, array_ID, x_pos, y_pos, z_pos, conn_type, periodic=[False,False,False]):

        self.array_ID = array_ID
        
        assert isinstance(conn_type, Connection), "conn_type must be of type Connection"

        # set number of compartments
        self.nx = len(x_pos)-1
        self.ny = len(y_pos)-1
        self.nz = len(z_pos)-1
        self.n_compartments = self.nx * self.ny * self.nz
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos
        self.box_len = [x_pos[-1]-x_pos[0],y_pos[-1]-y_pos[0],z_pos[-1]-z_pos[0]]
        self.periodic = periodic

        # initialize compartment dictionary
        self.compartments = {}
        for i in range(len(x_pos)-1):
            lx = x_pos[i+1]-x_pos[i]
            posx = (x_pos[i],x_pos[i+1])
            for j in range(len(y_pos)-1):
                ly = y_pos[j+1]-y_pos[j]
                posy = (y_pos[j],y_pos[j+1])
                for k in range(len(z_pos)-1):
                    lz = z_pos[k+1]-z_pos[k]
                    posz = (z_pos[k],z_pos[k+1])
                    sa = {'xy' : lx*ly, 'yz' : ly*lz, 'xz' : lx*lz}
                    self.compartments[(i,j,k)] = Compartment3D((i,j,k),
                                                             pos=[posx,posy,posz],
                                                             array_ID=self.array_ID,
                                                             surface_area=sa)
        # add connections
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    if i > 0:
                        # make connection to the left
                        self.compartments[(i,j,k)].connect(self.compartments[(i-1,j,k)],conn_type)
                    if i < self.nx-1:
                        # make connection to the right
                        self.compartments[(i,j,k)].connect(self.compartments[(i+1,j,k)],conn_type)
                    if j > 0:
                        # make connection to the bottom
                        self.compartments[(i,j,k)].connect(self.compartments[(i,j-1,k)],conn_type)
                    if j < self.ny-1:
                        # make connection to the top
                        self.compartments[(i,j,k)].connect(self.compartments[(i,j+1,k)],conn_type)
                    if k > 0:
                        # make connection to the bottom
                        self.compartments[(i,j,k)].connect(self.compartments[(i,j,k-1)],conn_type)
                    if k < self.nz-1:
                        # make connection to the top
                        self.compartments[(i,j,k)].connect(self.compartments[(i,j,k+1)],conn_type)

        if periodic[0]:
            for j in range(self.ny):
                for k in range(self.nz):
                    self.compartments[(0,j,k)].connect(self.compartments[(self.nx-1,j,k)],conn_type)
                    self.compartments[(self.nx-1,j,k)].connect(self.compartments[(0,j,k)],conn_type)

        if periodic[1]:
            for i in range(self.nx):
                for k in range(self.nz):
                    self.compartments[(i,0,k)].connect(self.compartments[(i,self.ny-1,k)],conn_type)
                    self.compartments[(i,self.ny-1,k)].connect(self.compartments[(i,0,k)],conn_type)

        if periodic[2]:
            for i in range(self.nx):
                for j in range(self.ny):
                    self.compartments[(i,j,0)].connect(self.compartments[(i,j,self.nz-1)],conn_type)
                    self.compartments[(i,j,self.nz-1)].connect(self.compartments[(i,j,0)],conn_type)

    def join3D(self,other_array,conn_type,append_side=None):
        """Joins a 3D compartment array to this
        one, making connections between compartments along
        one of the faces.

        conn_type is the Connection to be used to connect the 
        compartments between the two arrays.

        append_side must be one of [x-,x+,y-,y+,z-,z+], and denotes
        the face of "self" that will adjoin the other_array

        Note that if the other array is a CompartmentArray3D, then
        append_side also determines the face of the 
        other_array to be joined.  (e.g. if append_side = 'x-' then
        the x+ face of other_array will be joined to x-.
        """

        assert isinstance(other_array, CompartmentArray3D), "type of other_array must be CompartmentArray3D"

        assert append_side in ['x-','x+','y-','y+','z-','z+'], "append_side must be set to one of [x-,x+,y-,y+,z-,z+]"

        if append_side[0] == 'x':
            face1 = (self.ny,self.nz)
            face2 = (other_array.ny,other_array.nz)
        elif append_side[0] == 'y':
            face1 = (self.nx,self.nz)
            face2 = (other_array.nx,other_array.nz)
        elif append_side[0] == 'z':
            face1 = (self.nx,self.ny)
            face2 = (other_array.nx,other_array.ny)

        if face1 != face2:
            raise ValueError("Error! dimensions of arrays don't match along append dimensions: {0} {1}".format(face1,face2))

        assert isinstance(conn_type, Connection), "conn_type must be of type Connection"

        # add connections between compartments
        if append_side[0] == 'x':
            if append_side[1] == '-':
                xind1 = 0
                xind2 = other_array.nx-1
            else:
                xind1 = self.nx-1
                xind2 = 0
            for j in range(self.ny):
                for k in range(self.nz):
                    self.compartments[(xind1,j,k)].connect(other_array.compartments[(xind2,j,k)],conn_type)
                    other_array.compartments[(xind2,j,k)].connect(self.compartments[(xind1,j,k)],conn_type)
        elif append_side[0] == 'y':
            if append_side[1] == '-':
                yind1 = 0
                yind2 = other_array.ny-1
            else:
                yind1 = self.ny-1
                yind2 = 0
            for i in range(self.nx):
                for k in range(self.nz):
                    self.compartments[(i,yind1,k)].connect(other_array.compartments[(i,yind2,k)],conn_type)
                    other_array.compartments[(i,yind2,k)].connect(self.compartments[(i,yind1,k)],conn_type)
        elif append_side[0] == 'z':
            if append_side[1] == '-':
                zind1 = 0
                zind2 = other_array.nz-1
            else:
                zind1 = self.nz-1
                zind2 = 0
            for i in range(self.nx):
                for j in range(self.ny):
                    self.compartments[(i,j,zind1)].connect(other_array.compartments[(i,j,zind2)],conn_type)
                    other_array.compartments[(i,j,zind2)].connect(self.compartments[(i,j,zind1)],conn_type)

        
        
