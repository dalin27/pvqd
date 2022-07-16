# Projected - Variational Quantum Dynamics (p-VQD)

This repository contains the source code of the p-VQD algorithm presented in: 

Stephano Barison, Filippo Vicentini and Giuseppe Carleo, An efficient quantum algorithm for the time evolution
of parameterized circuits, 2021 (https://arxiv.org/abs/2101.04579).

## Content of the repository

- **ansatze.py** : creates the (fixed) quantum circuit ansatz; 
- **model.py** : prepares Pauli operators and the Hamiltonian. Here we implemented the Transverse Field Ising Model;
- **pVQD.py** : class of the pVQD algorithm, with methods to evaluate the overlap, the gradient, 
the observables and run;
- **example.py** : an example to use the pVQD method to simulate the Transverse Field Ising Model on an open chain of 3 
qubits;
- **plotter.py** : plot the results of the pVQD calculations and compare them to the exact classical simulation;
- **data** : folder that contains some pre-produced data to plot.

### Modified version of the code by Stefano Barison (https://github.com/StefanoBarison/p-VQD)
Here is a list of some things I changed:
- **model.py** (previously "pauli_function.py"): combined the functions "construct_total_circuit" and
"construct_total_circuit_local" in "pVQD.py" into one generalized function "tfim_hamiltonian" (to avoid redundancy);
- **pVQD.py** : modified the "projector_zero" and "projector_zero_local" functions to make them more compact; changed
the formatting of the "run" method (for instance by introducing the data log earlier);
- **example.py** : changed how to generate the observables to measure;

