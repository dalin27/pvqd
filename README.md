# Projected - Variational Quantum Dynamics (p-VQD)

This repository contains the source code of the p-VQD algorithm presented in: 

Stephano Barison, Filippo Vicentini and Giuseppe Carleo, An efficient quantum algorithm for the time evolution
of parameterized circuits, 2021 (https://arxiv.org/abs/2101.04579).

## Content of the repository

- **pVQD.py** : class of the pVQD algorithm, with methods to evaluate the overlap, the gradient, 
the observables and run.
- **model.py** : prepares Pauli operators and the Hamiltonian. Here we implemented the Transverse Field Ising Model.
- **ansatze.py** : creates the (fixed) quantum circuit ansatz. 
- **example.py** : an example to use the pVQD method to simulate the Transverse Field Ising Model on an open chain of 3 
qubits.
- **plotter.py** : plot the results of the pVQD calculations and compare them to the exact classical simulation.
- **data** : folder that contains some pre-produced data to plot.

### Modified version of the code by Stefano Barison (see the associated GitHub repository)
In particular, I updated some depreciated features of the code (ex: changed np.bool to bool), modified the 
"projector_zero" and "projector_zero_local" functions in "pVQD.py" to make them more compact, combined the functions 
"construct_total_circuit" and "construct_total_circuit_local" in "pVQD.py" into one generalized function (to avoid 
redundancy), changed the functions in "model.py" (previously called "pauli_function.py") to work directly with "I,X,Y,Z"
(and adapted accordingly the few lines of code in "example.py" defining the observables of interest), changed variable
names to make the code more consistent, and added comments throughout the code. 
