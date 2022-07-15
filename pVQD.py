# Import the relevant libraries
from qiskit.opflow import I, Z
from model import tens_prod
import numpy as np
from qiskit.opflow import CircuitSampler, StateFn
from qiskit.circuit import ParameterVector
from qiskit.opflow.evolutions import PauliTrotterEvolution
from qiskit.opflow.expectations import PauliExpectation, AerPauliExpectation
import json

# Define some utility functions


# Create a global projector |00...0><00...0|, which is required for post-selection of |0> qubit after applying
# C^\dagger(w+dw) e^{-i H \delta_t} C(w).
def projector_zero_global(num_qubits):
	return tens_prod([0.5*(I+Z)]*num_qubits)  # returns a PauliSumOp object


# Create a local version of the cost function proposed by Cerezo et al:
# https://www.nature.com/articles/s41467-021-21728-w.
def projector_zero_local(num_qubits):

	tot_prj = 0  # initialize the sum of all local projectors
	loc_prj_list = [I] * num_qubits  # initialize the local projector list

	for k in range(num_qubits):
		loc_prj_list[k] = 0.5*(I+Z)
		tot_prj += tens_prod(loc_prj_list)
		loc_prj_list[k] = I  # restore loc_prj_list to its initial form for the next iteration

	return tot_prj / num_qubits  # returns a PauliSumOp object


# Create a versor (to shift the parameters when using the parameter shift rule)
def ei(i,n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]


# This class simulates the time evolution of a quantum state by approximating it with a variational ansatz whose
# parameters are varied in order to follow the unitary evolution.
class pVQD:

	# Initialize variables and ParameterVector objects
	def __init__(self, hamiltonian, ansatz, ansatz_reps, init_params, init_shift, instance, shots):
		"""
		Args:
			hamiltonian (PauliSumOp class): Hamiltonian of interest for the evolution of the state
			ansatz (function with Qiskit QuantumCircuit class as output): fixed ansatz throughout the time evolution
			ansatz_reps (integer): how many times the ansatz is repeated (depth)
			init_params (list of real numbers): initialization of the variational parameters (w)
			init_shift (list of real numbers): initialization of the small shift of the variational parameters (dw)
			instance (Qiskit QuantumInstance class): quantum backend to execute circuits
			shots (integer): number of execution of Qiskit circuits
		"""

		self.hamiltonian = hamiltonian
		self.depth = ansatz_reps
		self.num_qubits = hamiltonian.num_qubits
		self.parameters = init_params
		self.num_params = len(init_params)
		self.shift = init_shift
		self.instance = instance
		self.shots = shots

		# Initialize a ParameterVector object to host all the variational parameters in the circuit ansatz
		self.params_circuit = ParameterVector('p', self.num_params)
		self.ansatz = ansatz(self.num_qubits, self.depth, self.params_circuit)

		# Initialize a ParameterVector object for the C^\dagger(w+dw) (left) and C(w) (right) unitaries.
		# Left meaning that it is applied last (if we draw the circuit it would appear on the far right)
		self.params_left = ParameterVector('l', self.num_params)
		self.params_right = ParameterVector('r', self.num_params)

		# Initialize a ParameterVector object for measuring observables
		self.params_obs = ParameterVector('θ', self.num_params)


	# Create the circuit that will be used to evaluate the overlap and its gradient
	def construct_circuit(self,time_step,projector_zero):
		"""
		Args:
			time_step (integer): small time step (for the Trotter step e^{-i H \delta_t})
			projector_zero (Qiskit PauliSumOp): either the global or local projector function defined above
		Returns:
			state_wfn (Qiskit ComposedOp): circuit state
		"""

		# Create the e^{-i H \delta_t} unitary
		step_h = time_step * self.hamiltonian
		trotter = PauliTrotterEvolution(reps=1) # "reps": how many Trotterization repetitions to make
		U_dt = trotter.convert(step_h.exp_i()).to_circuit()

		# Define the left and right circuits by assigning to the ansatz the associated parameters
		l_circ = self.ansatz.assign_parameters({self.params_circuit: self.params_left})
		r_circ = self.ansatz.assign_parameters({self.params_circuit: self.params_right})

		# Construct the circuit shown in Fig. 1 of the main Ref
		op_prj = StateFn(projector_zero(self.num_qubits), is_measurement=True) # get the projector (OperatorStateFn)
		# Implement C^\dagger(w+dw) e^{-i H \delta_t} C(w); Note C^\dagger = C^{-1}; returns a CircuitStateFn
		seq_gates = StateFn(r_circ + U_dt + l_circ.inverse())  # "+" deprecated; use "compose" instead
		state_wfn = op_prj @ seq_gates # circuit state (ComposedOp)

		return state_wfn


	# Calculate the overlap and its gradient using the parameter shift method (for gradient descent)
	def overlap_and_gradient_param_shift(self, state_wfn, parameters, shift, expectator, sampler):
		"""
		Args:
			state_wfn (Qiskit ComposedOp): circuit state
			parameters (numpy array of size N): current variational parameters (w)
			shift (numpy array of size N): current small shift in the variational parameters (dw)
			expectator (Expectation converter): expectations of quantum state circuits over Pauli observables
												(ex: PauliExpectation or AerPauliExpectation)
			sampler (CircuitSampler): uses a backend to convert any StateFn into an approximation of the state function
		Returns:
			overlap (numpy array with two entries): the overlap and its error
			gradient (numpy array of size Nx2): N gradients (one for each variational parameter) with its error
		"""

		# Initialize a list of dictionaries of the form "parameter: its value" for the overlap and the gradient.
		# The first dictionary corresponds to the parameters necessary to compute the overlap/loss L(dw) (Eq. 3).
		# This corresponds to only one circuit. Some "+" signs below join lists together (not performing addition).
		values_dict = [dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() +
							   (parameters + shift).tolist()))]

		# Append the results for the gradient (Eq. 6 in the main Ref, choosing s = \pi/2).
		# If there are N variational parameters, we need to measure the overlap on 2N different circuits:
		# one circuit for each of the shifted parameters dw+s*e_i (N of them) and dw-s*e_i (N of them).
		for i in range(self.num_params):
			values_dict.append(dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() + (parameters +
										shift + ei(i,self.num_params) * np.pi / 2.0).tolist())))
			values_dict.append(dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() + (parameters +
										shift - ei(i,self.num_params) * np.pi / 2.0).tolist())))

		results = []  # store the results of the evaluated circuits

		for values in values_dict:
			# Evaluate the circuits with the parameters assigned

			# todo: understand better how the measurements are performed

			sampled_op = sampler.convert(state_wfn,params=values)
			mean = sampled_op.eval().real
			est_err = 0  # error is 0 for a statevector-type simulator

			# If the backend is not a statevector-type simulator (not False = True), then do the following:
			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op).real
				est_err = np.sqrt(variance/self.shots)

			results.append([mean,est_err])

		overlap = np.zeros(2)
		gradient = np.zeros((self.num_params, 2))

		overlap[0], overlap[1] = results[0]

		for i in range(self.num_params):
			rplus = results[1+2*i]
			rminus = results[2+2*i]
			gradient[i,:] = (rplus[0]-rminus[0])/2.0, np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

		return overlap, gradient


	# Calculate the overlap and its gradient using the SPSA method
	def overlap_and_gradient_spsa(self, state_wfn, parameters, shift, expectator, sampler, count):

		# Define hyperparameters
		c = 0.1
		a = 0.16
		A = 1
		alpha = 0.602
		gamma = 0.101

		a_k = a/np.power(A+count,alpha)
		c_k = c/np.power(count,gamma)

		# Determine the random shift
		delta = np.random.binomial(1, 0.5, size=self.num_params)
		delta = np.where(delta == 0, -1, delta)
		delta = c_k*delta

		# First create the dictionary for overlap
		values_dict = [dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() +
							   (parameters + shift).tolist()))]

		# Then the values for the gradient
		values_dict.append(dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() +
								   (parameters + shift + delta).tolist())))
		values_dict.append(dict(zip(self.params_right[:] + self.params_left[:], parameters.tolist() +
								   (parameters + shift - delta).tolist())))

		results = []

		for values in values_dict:
			# Now evaluate the circuits with the parameters assigned
			sampled_op = sampler.convert(state_wfn,params=values)

			mean = sampled_op.eval()[0]
			# mean = np.power(np.absolute(mean),2)
			est_err = 0

			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op)[0].real
				est_err = np.sqrt(variance/self.shots)

			results.append([mean,est_err])

		overlap = np.zeros(2)
		gradient = np.zeros((self.num_params, 2))

		overlap[0], overlap[1] = results[0]

		rplus = results[1]
		rminus = results[2]

		for i in range(self.num_params):
			# G = (Ep - Em)/2Δ_i
			# var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
			gradient[i,:] = a_k*(rplus[0]-rminus[0])/(2.0*delta[i]),np.sqrt(rplus[1]**2+rminus[1]**2)/(2.0*delta[i])

		return overlap, gradient


	# This function calculates the expectation value of an operator w.r.t a circuit with fixed parameters
	def measure_op(self, obs_wfn, pauli_op, parameters, expectator, sampler):
		"""
		Args:
			obs_wfn (Qiskit QuantumCircuit object): state circuit
			pauli_op (PauliOp object): Pauli string of local operators
			parameters (numpy array of size N): current variational parameters (w)
			expectator (Expectation converter): expectations of quantum state circuits over Pauli observables
												(ex: PauliExpectation or AerPauliExpectation)
			sampler (CircuitSampler): uses a backend to convert any StateFn into an approximation of the state function
		Returns:
			res (list of two real numbers): result of the expectation value with its error
		"""

		# Prepare the operator and the parameters
		wfn = StateFn(obs_wfn)
		op = StateFn(pauli_op, is_measurement = True)
		values_obs = dict(zip(self.params_obs[:], parameters.tolist()))

		# Evaluate the circuit with the parameters assigned
		state_wfn = op @ wfn
		state_wfn = expectator.convert(state_wfn)
		sampled_op = sampler.convert(state_wfn, params=values_obs)
		mean_value = sampled_op.eval().real
		est_err = 0  # error is 0 for a statevector-type simulator

		# If the backend is not a statevector-type simulator (not False = True), then do the following:
		if (not self.instance.is_statevector):
			variance = expectator.compute_variance(sampled_op).real
			est_err = np.sqrt(variance/self.shots)

		res = [mean_value,est_err]

		return res

	# Run the pVQD algorithm
	def run(self, overlap_threshold, time_step, n_steps, obs_dict = {}, filename = 'data/output_algo.json',
		    max_iter = 100, opt = 'sgd', cost_fun = 'global', grad = 'param_shift'):
		"""
		Args:
			overlap_threshold (float): when the value is attained, the optimization of dw is deemed completed
			time_step (integer): small time step for the time evolution
			n_steps (integer): number of small time steps
			obs_dict (dictionary): observables to measure (of the form "name of op: PauliOp", where name of op is a str)
			filename (string): name of the file to output the results
			max_iter (integer): maximum number of iterations to optimize the small shift dw at a fixed time
			opt (string): optimizer (in principle 'sgd', 'adam', 'momentum', etc., but only 'sgd' is implemented here)
			cost_fun (string): cost function is either 'global' or 'local' to select the appropriate projector
			grad (string): method to compute the gradient ('param_shift' or 'spsa')
		Returns:
			file with all the data
		"""

		# Initialize the expectation converter and the circuit sampler
		if(self.instance.is_statevector and cost_fun == 'local'):
			expectation = AerPauliExpectation()
		else: 
			expectation = PauliExpectation()
		sampler = CircuitSampler(self.instance)

		# Prepare the circuit state to compute the overlap and its gradient
		if cost_fun == 'global':
			state_wfn = self.construct_circuit(time_step, projector_zero_global)
		elif cost_fun == 'local':
			state_wfn = self.construct_circuit(time_step, projector_zero_local)
		state_wfn = expectation.convert(state_wfn)

		# Prepare the circuit state to measure the observables of interest
		obs_wfn = self.ansatz.assign_parameters({self.params_circuit: self.params_obs})

		# List of all the timestamp of interest
		times = [i * time_step for i in range(n_steps)]

		# Create a data log to store all the results (to output in a file)
		data_log = {'num_iter': [],  # store the amount of iterations necessary to optimize dw for a fixed time step
					   'times': times,  # List of all the timestamp of interest
			 'init_fidelities': [],  # initial fidelity/overlap before optimizing dw (fixed time step)
		      'fin_fidelities': [],  # final fidelity/overlap after optimizing dw (fixed time step)
		        'init_err_fid': [],  # initial error on the fidelity/overlap before optimizing dw (fixed time step)
		         'fin_err_fid': [],  # final error on the fidelity/overlap after optimizing dw (fixed time step)
		              'params': [list(self.parameters)]  # updated variational parameters (fixed time step)
					}

		# If the dictionary of observables is not empty, do the following:
		if len(obs_dict) > 0:

			# Make preliminary measurements before starting the algorithm (time = 0)
			for (obs_name,obs_pauli) in obs_dict.items():
				first_measure = self.measure_op(obs_wfn, obs_pauli, self.parameters, expectation, sampler)
				data_log[str(obs_name)] = [first_measure[0]]
				data_log['err_'+str(obs_name)] = [first_measure[1]]

		print("\nStart of the algorithm")

		# Loop over the number of time steps (not up to the last element of times because of preliminary measurements)
		for i, time in enumerate(times[:-1]):

			print('\n================================== \n')
			print(f'Time slice: {i+1};    Time: {time}')
			# print("Shift before optimizing this step:",self.shift)
			# print("Initial parameters:", self.parameters)
			print('\n================================== \n')
			
			count = 0  # count the number of iterations to optimize the small shift dw at a fixed time
			overlap = [0.01,0]  # initialize the overlap value and its error at a fixed time

			# Iterate over the optimization steps of the small variational shift parameters dw for a fixed time
			while overlap[0] < overlap_threshold and count < max_iter:

				print("Shift optimizing step:", count + 1)

				count += 1  # add one iteration to the count (fixed time)

				# Measure energy and gradient
				if grad == 'param_shift':
					overlap, gradient = self.overlap_and_gradient_param_shift(state_wfn, self.parameters,
																			  self.shift, expectation, sampler)
				elif grad == 'spsa':
					overlap, gradient = self.overlap_and_gradient_spsa(state_wfn, self.parameters, self.shift,
																	   expectation, sampler, count)

				# Store the overlap in the log before optimizing dw
				if count == 1:
					data_log['init_fidelities'].append(overlap[0])
					data_log['init_err_fid'].append(overlap[1])

				# print('Overlap', overlap)
				# print('Gradient', gradient[:,0])

				if opt == 'sgd':
					self.shift = self.shift + gradient[:,0]

			print('\n---------------------------------- \n')
			# print("Shift after optimizing:", self.shift)
			# print("New parameters:", self.parameters + self.shift)
			print("New overlap: ", overlap[0])

			# Update the variational parameters for the next time step
			self.parameters = self.parameters + self.shift
				
			# Before the next time step, store the data in the log for that fixed time
			data_log['num_iter'].append(count)
			data_log['fin_fidelities'].append(overlap[0])
			data_log['fin_err_fid'].append(overlap[1])
			data_log['params'].append(list(self.parameters))
			if len(obs_dict) > 0:
				for (obs_name,obs_pauli) in obs_dict.items():
					run_measure = self.measure_op(obs_wfn, obs_pauli, self.parameters, expectation, sampler)
					data_log[str(obs_name)].append(run_measure[0])
					data_log['err_'+str(obs_name)].append(run_measure[1])

		# After going through all the time steps, print the following:
		print("\nEnd of the algorithm\n")
		print("Total measurements:", sum(data_log['num_iter']))
		print("Measure per step:", sum(data_log['num_iter'])/n_steps)

		# Write down the results in a file
		# todo: output data with a header and, generally, make the ouput file more readable
		json.dump(data_log, open(filename, 'w+'), indent = 6)
