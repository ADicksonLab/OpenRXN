"""
Connections govern the transport between compartments.
In some cases the transport can be described by a first-order
rate equation, e.g.:

d_n1/dt = -k12 * n1 + k21 * n2
d_n2/dt = -k21 * n2 + k12 * n1

where nX is the number of a given species in compartment X,
kAB is a rate constant governing transport between compartment 
A and B.

A note on units:
The k values in the above equation are in units of 1/sec,
and this is what is stored in the connection objects.
"""

from . import unit
import logging

class Connection(object):
    """Basis class for connections.  The argument species_rates
    is a dictionary, with Species IDs as keys and values hold 
    rate constant tuples (k12,k21)."""

class AnisotropicConnection(Connection):
    
    def __init__(self, species_rates,dim=3):
        """
        Directional (anisotropic) transport connection between two compartments. 

        Parameters
        ----------
        species_rates : dict
            Dictionary of Species, where the keys are Species IDs and the 
            values are tuples of transport rates (k_out,k_in).
            k_out corresponds to transport from compartment 1 to 2.
            k_in corresponds to transport from compartment 2 to 1.

                - A tuple (k_out, k_in) specifies directional transport rates.
                - A scalar quantity k is interpreted as symmetric transport
                and internally turn into tuples of transport rates (k, k).
                - All rates must be specified in units of 1/s.

        dim : int, optional
            Spatial dimension.

        """
        self.species_rates = species_rates
        self.dim = dim

        for s,r in self.species_rates.items():
            if not isinstance(r,tuple):
                logging.warning(f"Species {s}: one scalar rate provided. Assigning k_out == k_in.")
                r = (r,r)
                self.species_rates[s] = r

            elif len(r) != 2:
                raise ValueError(f"Species {s}: rate tuple must be length 2, got length {len(r)}.")

            self.species_rates[s][0].ito(1/unit.sec)
            self.species_rates[s][1].ito(1/unit.sec)

    @staticmethod
    def _flip_tuple(t):
        return (t[1],t[0])
            
    def reverse(self):
        rev_species_rates = {}
        for s,r in self.species_rates.items():
            rev_species_rates[s] = self._flip_tuple(r)

        return AnisotropicConnection(rev_species_rates)
            
class IsotropicConnection(Connection):

    def __init__(self, species_rates,dim=3):
        """
        Symmetric (isotropic) transport connection between two compartments. 

        Parameters
        ----------
        species_rates : dict
            Dictionary of Species, where the keys are Species IDs and the 
            value is a transport rate constant (k).

                - A scalar quantity k is interpreted as symmetric transport
                and internally turn into tuples of transport rates (k, k).
                - All rates must be convertible to units of 1/s.

        dim : int, optional
            Spatial dimension.

        """
        self.species_rates = species_rates
        self.dim = dim

        for s,r in self.species_rates.items():
            if not isinstance(k,tuple):
                k = (k,k)

            self.species_rates[s][0].ito(1/unit.sec)
            self.species_rates[s][1].ito(1/unit.sec)

class DivByVConnection(Connection):

    def __init__(self, species_rates,dim=3):
        """
        Volume-based transport connection between two compartments. 
        
        This connection stores transport coefficients proportional to compartment volume. 
        Rates should be specified in units of L^d/s, where L is length and d is the 
        compartment volume. 

        Parameters
        ----------
        species_rates : dict
            Dictionary of Species, where the keys are Species IDs and the 
            value is a transport rate constant (k).

        dim : int, optional
            Spatial dimension.

        """

        self.species_rates = species_rates
        self.dim = dim

        for s,k in self.species_rates.items():
            if not isinstance(k, tuple):
                k = (k,k)
                self.species_rates[s] = k
            self.species_rates[s][0].ito(unit.nm**self.dim/unit.sec)
            self.species_rates[s][1].ito(unit.nm**self.dim/unit.sec)
            
class FicksConnection(Connection):

    def __init__(self, species_d_constants, surface_area=None, ic_distance=None, dim=3):
        """FicksConnection types use diffusion constants for each
        Species, together with the widths and adjoining surface area
        of the compartments, to determine rate constants for transport.

        F_net = D * A * DeltaC / DeltaX = F21 - F12

        where F_net is the net diffusive flux from compartment 2 to 
        compartment 1 (in molecules per unit time), 
        D is the diffusion constant (for a given species),
        A is the adjoining surface area between the compartments,
        DeltaX is the distance between the compartment centers.

        This is separated into two first order diffusion reactions:
        
        F21 = D * A * C2/DeltaX = n2*kV21 / V2
        F12 = D * A * C1/DeltaX = n1*kV12 / V1

        where kV12 = D*A/DeltaX and kV21 = D*A/DeltaX are equal to the
        compartment volumes times the first order rate constants of 
        diffusion between compartment 1 and compartment 2 (in 1/s).  
        kV12 and kV21 will be used to construct a DivByVConnection 
        object when the FicksConnection is resolved.

        surface_area is the surface area of the connecting interface between 
        the compartments (optional)

        ic_distance is the distance between the centers of the compartments
        (optional)

        If either surface_area or ic_distance is left undefined, they will be
        automatically calculated using compartment positions.

        Parameters
        ----------
        species_d_constants : dict
            Dictionary of Species, where the keys are Species IDs and the 
            value is a diffusion constant (D). (units convertible to L^2 / s).
        surface_area : Quantity, optional
            Interface area A between compartments (units convertible to L^(dim-1)).
        ic_distance : Quantity, optional
            Center-to-center distance DeltaX (units convertible to L).
        dim : int, optional
            Spatial dimension (default 3).

        """

        self.species_d_constants = species_d_constants
        self.surface_area = surface_area
        self.ic_distance = ic_distance
        self.dim = dim

    def resolve(self):
        """This returns an IsotropicConnection that does not 
        require any information about the Species, or the arrays"""

        if self.surface_area is None or self.ic_distance is None:
            raise ValueError("Error! This connection is not ready to be resolved.")

        rates = {}
        for s,d in self.species_d_constants.items():
            kV = d*self.surface_area/self.ic_distance
            kV.ito(unit.nm**self.dim/unit.sec)
            rates[s] = kV
            
        return DivByVConnection(rates,self.dim)

class ResConnection(Connection):

    def __init__(self, species_d_constants, surface_area=None, ic_distance=None, dim=3, face=None):
        """
        ResConnections are special connections from a compartment to a 
        reservoir.  They are resolved similarly to a FicksConnection.

        Parameters
        ----------
        species_d_constants : dict
            Dictionary of Species, where the keys are Species IDs and the 
            value is a diffusion constant (D). (units convertible to L^2/s).
        surface_area : Quantity, optional
            Reservoir interface area A (units convertible to L^(dim-1)).
        ic_distance : Quantity, optional
            Distance DeltaX from compartment center to reservoir boundary (units convertible to L).
        dim : int, optional
            Spatial dimension.
        face : str({'x','y','z'}) or None, optional
            Axis normal to the reservoir interface.

        """
        
        self.species_d_constants = species_d_constants
        self.surface_area = surface_area
        self.ic_distance = ic_distance
        self.dim = dim
        self.face = face

    def resolve(self):
        """This returns an AnisotropicConnection that does not 
        require any information about the Species, or the arrays"""

        if self.surface_area is None or self.ic_distance is None:
            raise ValueError("Error! This connection is not ready to be resolved.")

        rates = {}
        for s,d in self.species_d_constants.items():
            kV = d*self.surface_area/self.ic_distance
            kV.ito(unit.nm**self.dim/unit.sec)
            rates[s] = kV
            
        return DivByVConnection(rates,self.dim)
