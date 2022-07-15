# Import the relevant scripts and libraries
from ansatz import *
from model import *
from pVQD import *
from qiskit import Aer
from qiskit.utils import QuantumInstance

# todo: generate strings of operators in the Hamiltonian using itertools product or combinations (for operator pool)

# Initialize system parameters for Ising
num_qubits = 3  # number of qubits/spins
coup = 0.25  # coupling constant
field = 1.0  # field strength
dt = 0.05  # small time step
n_steps = 40  # number of time steps

# Algorithm parameters
overlap_treshold = 0.99999  # threshold for the overlap (the overlap refers to Eq. 1 in the main Ref)
depth = 3  # depth of the circuit
dw = 0.01  # initial shift

# Initialize the variational parameters and the associated circuit.
# 0's as initial parameters imply that the circuit ansatz is initially the identity operator (C(0)=1)
init_params = np.zeros((depth+1) * num_qubits + depth * (num_qubits-1))
wfn = hweff_ansatz(num_qubits, depth, init_params)

# Generate the Hamiltonian
H = tfim_hamiltonian(num_qubits, coup, field)

print(f'Initial circuit:\n{wfn}')
print(f'Hamiltonian:\n{H}')

# Set the initial shift
shift = np.array([dw]*len(init_params))
print(f'Initial shift:\n{shift}')

# Qiskit Aer currently includes three high performance simulator backends. The two of interest:
# QasmSimulator: Allows ideal and noisy multi-shot execution of qiskit circuits and returns counts or memory
# StatevectorSimulator: Allows ideal single-shot execution of qiskit circuits and returns the final statevector of
# the simulator after application
shots = 800
backend = Aer.get_backend('statevector_simulator')  # can replace to 'qasm_simulator'
instance = QuantumInstance(backend=backend, shots=shots)

obs = {}  # Initialize a dictionary to store the observables to measure
loc_op = [I] * num_qubits  # Initialize a list where each entry corresponds to a local operator acting on qubit i

# Magnetization operators
for i, k in enumerate(reversed(range(num_qubits))):
	loc_op[k] = Z; obs['Sz_'+str(i)] = tens_prod(loc_op)
	loc_op[k] = X; obs['Sx_'+str(i)] = tens_prod(loc_op)
	loc_op[k] = Y; obs['Sy_'+str(i)] = tens_prod(loc_op)
	loc_op[k] = I  # restore loc_op to its initial form for the next iteration

print('Observables:')
for (name_op, pauli_str) in obs.items():
	print(f'{name_op} --> {pauli_str}')

# Define the necessary parameters to run the algorithm
opt = 'sgd'  # choose gradient optimizer: 'sgd', 'adam', 'momentum' (only 'sgd' implemented here)
cost = 'local'  # choose which type of cost function to use: 'global', 'local'
grad = 'param_shift'  # choose how to estimate the gradient on hardware: 'param_shift', 'spsa'

# Instantiate the pVQD class
algo = pVQD(hamiltonian   = H,
			ansatz        = hweff_ansatz,
			ansatz_reps   = depth,
			init_params   = init_params,
			init_shift    = shift,
			instance      = instance,
			shots         = shots)

# Run the algorithm
algo.run(overlap_treshold,
		 dt,
		 n_steps,
	     obs_dict      = obs,
	     filename      = 'data/tfim_nspins={}_j={}_h={}_statevec.json'.format(num_qubits,coup,field),
	     max_iter      = 50,
	     opt           = opt,
	     cost_fun      = cost,
	     grad          = grad)