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

import logging

class Species(object):
    """Species can be either reactants or products and can
    interconvert through Reactions. Each species object must 
    be given a unique name."""

    def __init__(self, ID, **kwargs):
        
        self.ID = ID
        for key, value in kwargs.items():
            self.key = value

class Reaction(object):
    """A reaction describes how reactants and products are 
    related
    """

    def __init__(self, ID, reactants, products, stoich_r, stoich_p, kf=0, kr=0):
        self.ID = ID

        for r in reactants:
            assert isinstance(r,Species), "Error! Reactants must be Species objects"
        
        self.reactants = reactants
        assert len(reactants) == len(stoich_r), "Error! Stoichiometry list for reactants must have same length as reactants list"
        self.stoich_r = stoich_r

        for p in products:
            assert isinstance(p,Species), "Error! Products must be Species objects"
        
        self.products = products
        assert len(products) == len(stoich_p), "Error! Stoichiometry list for products must have same length as products list"
        self.stoich_p = stoich_p

        if kf < 0 or kr < 0:
            raise ValueError("Error!  Reaction rate cannot be negative")
            
        if kf == 0 and kr == 0:
            logging.warn("Warning: Both forward and reverse rates are set to zero")

        self.kf = kf
        self.kr = kr

        self.reactant_IDs = [s.ID for s in self.reactants]
        self.product_IDs = [s.ID for s in self.products]

        # todo: assure that the units on the rates are correct

    def display(self):
        """Returns a print string summarizing the reaction."""
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

        return(to_print + rate_str)

        
