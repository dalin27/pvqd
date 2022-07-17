# Import the relevant libraries
from qiskit.opflow import I, X, Y, Z
from functools import reduce

# Compute the tensor product of a list of PauliOp
def tens_prod(op_list):
	return reduce(lambda x, y: x ^ y, op_list)  # "x ^ y" is equivalent to "x.tensor(y)"


# Compute sigma_i = I^{i-1} \otimes sigma \otimes I^{N-i}, where sigma is a Pauli operator
def one_site_op(loc_op, site, pauli_op):
	"""
	Args:
		loc_op (list of size num_qubits): list of (local) identity operator I
		site (integer): index where the operator pauli_op acts non-trivially; site=0,...,n_spins-1
		pauli_op (PauliOp): either I,X,Y,Z
	Returns:
		PauliOp object acting on all n_spins
	"""
	loc_op[-site-1] = pauli_op  # "-site" because PauliOp string read from right to left; "-1" because site starts at 0
	tot_op = tens_prod(loc_op)
	loc_op[-site-1] = I  # restore loc_op to its initial form
	return tot_op


# Compute sigma_i * sigma_{i+1}= I^{i-1} \otimes sigma \otimes sigma \otimes I^{N-i-1}
def nearest_neigh_op(loc_op, site, pauli_op):
	"""
	Args:
		loc_op (list of size num_qubits): list of (local) identity operator I
		site (integer): index where the operator pauli_op acts non-trivially; site=0,...,n_spins-1
		pauli_op (PauliOp): either I,X,Y,Z
	Returns:
		PauliOp object acting on all n_spins
	"""
	if site != len(loc_op)-1:  # if "site" is not the last spin index [len(loc_op) = n_spins-1], do the following:
		loc_op[-site-1] = loc_op[-site-2] = pauli_op  # "-site" because PauliOp string read from right to left
		tot_op = tens_prod(loc_op)
		loc_op[-site-1] = loc_op[-site-2] = I  # restore loc_op to its initial form

	elif site == len(loc_op)-1:  # if "site" is the last spin index <==> if the boundary conditions are periodic
		loc_op[0] = loc_op[-1] = pauli_op  # add one pauli_op at both ends of the list
		tot_op = tens_prod(loc_op)
		loc_op[0] = loc_op[-1] = I  # restore loc_op to its initial form

	return tot_op


# Implement the Transverse Field Ising Model (TFIM) with either periodic or open boundary conditions
def tfim_hamiltonian(n_spins, model_params, bound_cond = 'open'):
	"""
	Args:
		n_spins (integer): number of num_qubits/qubits
		model_params (list of 2 floats): coupling parameter and field strength in a list
		bound_cond (string): boundary condition; either 'periodic' or 'open'
	Returns:
		h (PauliSumOp object): Hamiltonian of the transverse field Ising model
	"""

	coup, field = model_params  # unpack the model parameters
	h = 0  # initialize the Hamiltonian
	loc_op = [I] * n_spins # initialize a list where each entry corresponds to a local operator acting on qubit i

	# X terms
	for k in range(n_spins):
		h += field * one_site_op(loc_op, k, X)

	# ZZ terms
	n_max = n_spins-1 if bound_cond == 'open' else n_spins  # if bound_cond == 'periodic', set n_max = n_spins
	for k in range(n_max):
		h += coup * nearest_neigh_op(loc_op, k, Z)

	return h


# Implement the Heisenberg XYZ model with either periodic or open boundary conditions
def xyz_hamiltonian(n_spins, model_params, bound_cond='open'):
	"""
	Args:
		n_spins (integer): number of num_qubits/qubits
		model_params (list of 4 floats): 3 coupling parameters and the field strength
		bound_cond (string): boundary condition; either 'periodic' or 'open'

	Returns:
		h (PauliSumOp object): Hamiltonian of the Heisenberg XYZ model
	"""

	jx, jy, jz, field = model_params  # unpack the model parameters
	h = 0  # initialize the Hamiltonian
	loc_op = [I] * n_spins # initialize a list where each entry corresponds to a local operator acting on qubit i

	# Z terms
	for k in range(n_spins):
		h += field * one_site_op(loc_op, k, Z)

	# Interaction terms
	n_max = n_spins-1 if bound_cond == 'open' else n_spins  # if bound_cond == 'periodic', set n_max = n_spins
	for k in range(n_max):
		h += jx * nearest_neigh_op(loc_op, k, X)  # XX terms
		h += jy * nearest_neigh_op(loc_op, k, Y)  # YY terms
		h += jz * nearest_neigh_op(loc_op, k, Z)  # ZZ terms
		
	return h

