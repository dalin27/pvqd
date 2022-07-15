# Import the relevant libraries
from qiskit import QuantumCircuit

# Fixed ansatz throughout the time evolution (Eqs. 14-15 in the main Ref)
def hweff_ansatz(n_spins,depth,p):
	"""
	Args:
		n_spins (integer): number of qubits/num_qubits (N in the main Ref)
		depth (integer): depth of the circuit (d in the main Ref)
		p (vector of real numbers): variational parameters (w in the main Ref)
	Returns:
		circuit: Qiskit QuantumCircuit class
	"""

	circuit = QuantumCircuit(n_spins)
	count = 0

	# As in Eqs. 14-15, we take the product of the blocks from l=1 to d
	for l in range(depth):

		if(l % 2 == 0):
			# Rx - Rzz block
			for i in range(n_spins):
				circuit.rx(p[count],i)
				count += 1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count += 1

			circuit.barrier()

		if(l % 2 == 1):
			for i in range(n_spins):
				circuit.ry(p[count],i)
				count += 1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count += 1

			circuit.barrier()

	# Final block to close the ansatz
	if (depth % 2 == 1):
		for i in range(n_spins):
				circuit.ry(p[count],i)
				count += 1
	if (depth % 2 == 0):
		for i in range(n_spins):
				circuit.rx(p[count],i)
				count += 1

	return circuit

















