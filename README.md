# OpenRXN
A free, open-source tool for modeling chemical reaction networks in Python.  Offers a single platform for pharmacokinetics, pharmacodynamics, and reaction-diffusion equations. 

## Features
* definition of compartments, reactions and species as Python objects
* capable of stochastic (e.g. Gillespie) or deterministic modeling (e.g. ODEs)
* uses Pint for units throughout
* easy definition of 1D, 2D or 3D compartment arrays
* can visualize compartment connectivity as graphs (NetworkX)

## Dependencies
* NetworkX
* Numpy
* Scipy
* Pint
* pandas
* matplotlib

## Installation with pip

First, clone the source from github:

```
git clone https://github.com/ADicksonLab/OpenRXN.git
```

Change into that directory and install using pip:

```
cd OpenRXN
pip install .
```

Check the install by running one of the examples.  E.g.:
```
python examples/1D_diffusion.py
```

