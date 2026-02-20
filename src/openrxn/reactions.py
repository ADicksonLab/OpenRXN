"""Reaction objects describe the relationship between a set
of reactants and a set of products.  Both the reactants and 
products are Species objects.

The combustion reaction  CH4  +  2O2 -->  CO2  + 2H2O  
is represented as follows:

from openrxn.reactions import Reaction, Species
from openrxn import unit

m = Species('methane')
o2 = Species('molecular oxygen')
co2 = Species('carbon dioxide')
w = Species('water')
kf = 1000/unit.sec
kf = 0/unit.sec

r1 = Reaction('combustion',[m,o2],[co2,w],[1,2],[1,2],kf=kf,kr=kr)
"""

from openrxn import unit
import logging

class Species(object):
    """
    Species object models an individual chemical entity (reactants or product) that can participate and interconvert
    through Reactions. Each species object must be given a unique name ('ID').

    Parameters
    ----------
    ID : str
        Unique identifier for the species.
    **kwargs : dict, optional
        Arbitrary key-value pairs defining additional species properties.
    """

    def __init__(self, ID, **kwargs):
        self.ID = ID

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return f"Species({self.ID})"

class Reaction(object):
    """
    A reaction describes how reactants and products are related.

    Parameters
    ----------
    ID: str
        Unique identifier for the species.
    reactants: set of Species
        Reactant Species objects participate in the reaction
    products: set of Species
        Products Species objects participate in the reaction
    stoich_r: int
        Stoichiometric coefficients for reactants. Must be positive integers
    stoich_p: int
        Stoichiometric coefficients for products. Must be positive integers
    kf: float, optional
        Forward rate constant. 
    kr: float, optional
        Reverse rate constant. 
    """

    def __init__(self, ID, reactants, products, stoich_r, stoich_p, kf=0, kr=0):
        self.ID = ID

        if len(reactants) != len(stoich_r):
            raise ValueError("Stoichiometry list for reactants must have same length as reactants list")
        if len(products) != len(stoich_p):
            raise ValueError("Stoichiometry list for products must have same length as products list")

        if not all(isinstance(r, Species) for r in reactants):
            raise TypeError("Reactants must be Species objects")
        if not all(isinstance(p, Species) for p in products):
            raise TypeError("Products must be Species objects")

        if kf < 0 or kr < 0:
            raise ValueError("Error!  Reaction rate cannot be negative")

        if kf == 0 and kr == 0:
            logging.warning("Both forward and reverse rates are set to zero")

        self.reactants = reactants
        self.stoich_r = stoich_r
        self.products = products
        self.stoich_p = stoich_p
        self.kf = kf
        self.kr = kr 

        self.reactant_IDs = [s.ID for s in self.reactants]
        self.product_IDs = [s.ID for s in self.products]

    def _expected_unit_from_order(self, stoich):

        order = float(sum(stoich))
        power = 1.0 - order

        molar_factor = 1 if abs(power) < 1e-12 else (unit.molar ** power)
        return molar_factor / unit.second
    
    def _validate_rate(self, rate, expected_unit):
        """
        Assure that the units on the rates are correct.
        """

        if isinstance(rate, (int,float)):
            q = rate * expected_unit
            logging.warning(f"Reaction rate provided without units. Assigning {expected_unit}. ")
        else:
            q = unit.Quantity(rate)

        if q.magnitude < 0: 
            raise ValueError(f"Reaction rate cannot be negative!")
            
        if q.dimensionality != expected_unit.dimensionality:
            raise ValueError(f"Reaction rate has wrong units: {q.units}. "
                            f"Expected unit: {expected_unit}")
            
        return q.to(expected_unit)

    def __repr__(self):
        rxn_list = self.display()
        return f"Reaction({rxn_list})"

    def display(self):
        """
        Returns a string summarizing the reaction.
        """

        to_print = ""
        rate_str = ""
        
        for i,r in enumerate(self.reactants):
            if i > 0:
                to_print += "+ "
            if self.stoich_r[i] > 1:
                to_print += "{1} {0} ".format(r.ID,self.stoich_r[i])
            else:
                to_print += "{0} ".format(r.ID)

        if self.kr > 0:
            to_print += "<"
            rate_str += " // kr = {0}".format(self.kr)

        to_print += "---"

        if self.kf > 0:
            to_print += "> "
            rate_str += " // kf = {0}".format(self.kf)

        for i,p in enumerate(self.products):
            if i > 0:
                to_print += "+ "
            if self.stoich_p[i] > 1:
                to_print += "{1} {0} ".format(p.ID,self.stoich_p[i])
            else:
                to_print += "{0} ".format(p.ID)

        return (to_print + rate_str)


