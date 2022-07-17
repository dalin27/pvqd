# Import the relevant scripts and libraries
from ansatz import *
from model import *
from pVQD import *
from qiskit import Aer
from qiskit.utils import QuantumInstance

#----------------------------------------------------------------------------------------------------#
# # Choose the parameters for the Transverse Field Ising Model
# model = 'tfim'
# num_qubits = 3  # number of qubits/spins
# depth = 3  # depth of the circuit ansatz
# model_params = [0.25, 1.0]  # [j,h]; j: coupling constant, h: field strength

# Choose the parameters for the Heisenberg XYZ model
model = 'xyz'
num_qubits = 3  # number of qubits/spins
depth = 3  # depth of the circuit ansatz
model_params = [0.1, 0.1, 0.1, 1.0]  # [jx,jy,jz,h]; j: coupling constants, h: field strength

# Choose the algorithm parameters
dt = 0.05  # small time step
n_steps = 40  # number of time steps
overlap_treshold = 0.99999  # threshold for the overlap (the overlap refers to Eq. 1 in the main Ref)
dw = 0.01  # initial shift of the variational parameters
shots = 1  # number of times to execute Qiskit circuit (1 ==> statevector simulator, >1 ==> qasm simulator)
opt = 'sgd'  # choose gradient optimizer: 'sgd', 'adam', 'momentum' (only 'sgd' implemented here)
cost = 'local'  # choose which type of cost function to use: 'global', 'local'
grad = 'param_shift'  # choose how to estimate the gradient on hardware: 'param_shift', 'spsa'

# Create a dictionary to store all the parameters chosen by the user (to output in the result file)
def param_dictio(*args):
    return dict(((param_name, eval(param_name)) for param_name in args))
input_dictio = param_dictio('model','num_qubits','depth','model_params','dt','n_steps','overlap_treshold',
							'dw','shots','opt','cost','grad')
#----------------------------------------------------------------------------------------------------#


if model == 'tfim':

	# Generate the Hamiltonian
	H = tfim_hamiltonian(num_qubits, model_params)

	# Initialize the variational parameters and the circuit ansatz.
	# 0's as initial parameters imply that the circuit ansatz is initially the identity operator (C(0)=1)
	init_params = np.zeros((depth + 1) * num_qubits + depth * (num_qubits - 1))
	ansatz = tfim_ansatz
	wfn = ansatz(num_qubits, depth, init_params)

elif model == 'xyz':

	# Generate the Hamiltonian
	H = xyz_hamiltonian(num_qubits, model_params)

	# Initialize the variational parameters and the circuit ansatz.
	init_params = np.zeros((depth + 1) * num_qubits + 3 * depth * (num_qubits - 1))
	ansatz = xyz_ansatz
	wfn = ansatz(num_qubits, depth, init_params)

else:
	print('This Hamiltonian is not yet implemented.')

print(f'Initial circuit:\n{wfn}')
print(f'Hamiltonian:\n{H}')

# Set the initial shift
shift = np.array([dw]*len(init_params))
print(f'Initial shift:\n{shift}')

# Choose backend. Qiskit Aer currently includes three high performance simulator backends. The two of interest:
# QasmSimulator: Allows ideal and noisy multi-shot execution of Qiskit circuits and returns counts or memory
# StatevectorSimulator: Allows ideal single-shot execution of Qiskit circuits and returns the final statevector of
# the simulator after application.
if shots > 1:
	backend = Aer.get_backend('qasm_simulator')
	out_filename = 'data/test_{}_nspins={}_params='.format(model,num_qubits) \
				   + ('_{}' * len(model_params)).format(*model_params) + '_shots={}.json'.format(shots)
elif shots == 1:
	backend = Aer.get_backend('statevector_simulator')
	out_filename = 'data/test_{}_nspins={}_params='.format(model,num_qubits) \
				   + ('_{}' * len(model_params)).format(*model_params) + '_statevec.json'
else:
	print('This backend is not available.')
instance = QuantumInstance(backend=backend, shots=shots)

# Preparing observables
obs = {}  # Initialize a dictionary to store the observables to measure
loc_op = [I] * num_qubits  # Initialize a list where each entry corresponds to a local operator acting on qubit i

# Local magnetization operators
for k in range(num_qubits):
	obs['Sz_'+str(k)] = one_site_op(loc_op, k, Z)
	obs['Sx_'+str(k)] = one_site_op(loc_op, k, X)
	obs['Sy_'+str(k)] = one_site_op(loc_op, k, Y)

print('Observables:')
for (name_op, pauli_str) in obs.items():
	print(f'{name_op} --> {pauli_str}')

# Instantiate the pVQD class
algo = pVQD(hamiltonian   = H,
			ansatz        = ansatz,
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
	     filename      = out_filename,
	     max_iter      = 50,
	     opt           = opt,
	     cost_fun      = cost,
	     grad          = grad,
		 input_dictio  = input_dictio)