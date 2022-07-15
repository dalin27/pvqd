# Import the relevant libraries
from qiskit.opflow import I, X, Y, Z
from functools import reduce

def tens_prod(op_list):
	# take the tensor product of a list of PauliOp
	# "x ^ y" is equivalent to "x.tensor(y)"
	return reduce(lambda x, y: x ^ y, op_list)


def tfim_hamiltonian(n_spins, coup, field, bound_cond='open'):
	"""
	Args:
		n_spins (integer): number of num_qubits/qubits
		coup (float): coupling parameter
		field (float): field strength
		bound_cond (string): boundary condition; either 'periodic' or 'open' (or any other string that is not 'periodic')

	Returns:
		h (PauliSumOp object): Hamiltonian of the transverse field Ising model
	"""

	h = 0  # initialize the Hamiltonian

	# Initialize a list where each entry corresponds to a local operator acting on qubit i
	loc_op = [I] * n_spins

	# ZZ terms
	for k in range(n_spins - 1):
		loc_op[k] = loc_op[k + 1] = Z
		h += coup * tens_prod(loc_op)
		loc_op[k] = loc_op[k + 1] = I  # restore loc_op to its initial form for the next iteration

	# ZZ term for periodic boundary conditions
	if bound_cond == 'periodic':
		loc_op[0] = loc_op[-1] = Z
		h += coup * tens_prod(loc_op)
		loc_op[0] = loc_op[-1] = I  # restore loc_op to its initial form for the next iteration

	# X terms
	for k in range(n_spins):
		loc_op[k] = X
		h += field * tens_prod(loc_op)
		loc_op[k] = I  # restore loc_op to its initial form for the next iteration

	return h


# todo: implement the Heisenberg model (what would be a good ansatz?)