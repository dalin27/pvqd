# Import the relevant libraries
import json
import matplotlib.pyplot as plt


# Plot the results obtained with the pVQD algorithm and compare them to the exact classical simulations.
# In particular, we plot the expectation value of \sigma_x and \sigma_z measured on the first spin (wlog).

exact = json.load(open('data/exact_result_J0.25_B1.dat'))
data = json.load(open('data/tfim_nspins=3_params=_0.25_1.0_statevec.json'))

spin = '1'  # spin of interest

# Instantiate the plot
fig, ax = plt.subplots(2, 1, sharex=True)

# Plot the Sx measurements
ax[0].plot(exact['times'][:40],exact['Sx'][:40], label ="Exact", linestyle='dashed', linewidth=1.2, color='black')
ax[0].errorbar(data['times'][:40], data['Sx_'+str(spin)][:40], yerr=data['err_Sx_'+str(spin)][:40], label="pVQD",
               marker='o', linestyle='', elinewidth=1, color='C0', capsize=2, markersize=1.5)
ax[0].set(ylabel = r"$\langle\sigma_{x}\rangle_{1}$")
ax[0].set_ylim(ymax=1.1, ymin=-1.1)

# Plot the Sz measurements
ax[1].plot(exact['times'][:40], exact['Sz'][:40], label ="Exact", linestyle='dashed', linewidth=1.2, color='black')
ax[1].errorbar(data['times'][:40], data['Sz_'+str(spin)][:40], yerr=data['err_Sz_'+str(spin)][:40], label="pVQD",
               marker='o', linestyle='', elinewidth=1, color='C0', capsize=2, markersize=1.5)
ax[1].set(ylabel=r"$\langle\sigma_{z}\rangle_{1}$", xlabel=r'$t$')
ax[1].set_ylim(ymax=1.1, ymin=-1.1)

# Legend above the plots
lines, labels = ax[0].get_legend_handles_labels()
ax[0].legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=2, fancybox=True, shadow=False)

plt.show()
