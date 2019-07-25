"""A collection of Models that can be used as templates.
These mostly share the same features but differ in their 
initialization of compartments."""

from openrxn.model import Model
from openrxn.compartments import Compartment

class OneCompartmentModel(Model):
    """Uses a single compartment called main.  No arguments
    required."""
    def __init__(self):
        self.compartments = {'main' : Compartment('main')}

