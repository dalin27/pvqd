# Import the relevant libraries
from qiskit import QuantumCircuit

# Fixed ansatz throughout the time evolution for the TFIM (Eqs. 14-15 in the main Ref)
def tfim_ansatz(n_spins, depth, params):
	"""
	Args:
		n_spins (integer): number of qubits/num_qubits (N in the main Ref)
		depth (integer): depth of the circuit (d in the main Ref)
		params (vector of real numbers): variational parameters (w in the main Ref)
	Returns:
		circuit: Qiskit QuantumCircuit class
	"""

	circuit = QuantumCircuit(n_spins)
	count = 0

	# As in Eqs. 14-15, we take the product of the blocks from l=1 to d
	for l in range(depth):

		# For the even l-th block, apply rx and rzz gates
		if(l % 2 == 0):
			# Rx - Rzz block
			for i in range(n_spins):
				circuit.rx(params[count],i)
				count += 1
			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(params[count],i,i+1)
				count += 1
			circuit.barrier()

		# For the odd l-th block, apply ry and rzz gates
		if(l % 2 == 1):
			for i in range(n_spins):
				circuit.ry(params[count],i)
				count += 1
			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(params[count],i,i+1)
				count += 1
			circuit.barrier()

	# Final block to close the ansatz
	if (depth % 2 == 1):
		for i in range(n_spins):
				circuit.ry(params[count],i)
				count += 1
	if (depth % 2 == 0):
		for i in range(n_spins):
				circuit.rx(params[count],i)
				count += 1

	return circuit


# Fixed ansatz throughout the time evolution for the Heisenberg XYZ model
def xyz_ansatz(n_spins, depth, params):
	"""
	Args:
		n_spins (integer): number of qubits/num_qubits (N in the main Ref)
		depth (integer): depth of the circuit (d in the main Ref)
		params (vector of real numbers): variational parameters (w in the main Ref)
	Returns:
		circuit: Qiskit QuantumCircuit class
	"""

	circuit = QuantumCircuit(n_spins)
	count = 0

	for l in range(depth):

		# Rz, Rxx, Ryy, Rzz block
		for i in range(n_spins):
			circuit.rz(params[count],i)
			count += 1
		circuit.barrier()

		for i in range(n_spins-1):
			circuit.rxx(params[count],i,i+1)
			count += 1
		circuit.barrier()

		for i in range(n_spins-1):
			circuit.ryy(params[count],i,i+1)
			count += 1
		circuit.barrier()

		for i in range(n_spins-1):
			circuit.rzz(params[count],i,i+1)
			count += 1
		circuit.barrier()

	# Final block to close the ansatz
	for i in range(n_spins):
			circuit.rz(params[count],i)
			count += 1

	return circuit













