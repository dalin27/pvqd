# Projected Variational Quantum Dynamics (pVQD)

This repository is a modified version of the code by Stefano Barison (https://github.com/StefanoBarison/p-VQD), which
is based on the following article:

Stephano Barison, Filippo Vicentini and Giuseppe Carleo, An efficient quantum algorithm for the time evolution
of parameterized circuits, 2021 (https://arxiv.org/abs/2101.04579).

The skeleton of Stefano's code is the same, but some functions were added, modified and also removed. The layout was
also changed and comments have been added.

## Content of the repository

- **ansatze.py** : creates the (fixed) quantum circuit ansatz; 
- **model.py** : prepares Pauli operators and the Hamiltonian. Here we implemented the Transverse Field Ising Model 
(TFIM) and the Heisenberg XYZ model;
- **pVQD.py** : class of the pVQD algorithm, with methods to evaluate the overlap, the gradient, 
the observables and run;
- **example.py** : an example of code to run the pVQD method to simulate either the TFIM or the XYZ model;
- **plotter.py** : plot the results of the pVQD calculations and compare them to the exact classical simulation;
- **data** : folder that contains some pre-produced data to plot.


