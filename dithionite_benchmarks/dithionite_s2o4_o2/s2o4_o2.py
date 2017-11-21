# Id: abiotic.py, Mon 28 Aug 2017 05:03:53 PM MDT #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for abiotic reactions.
#------------------------------------------------------------------------------
import sys
import os
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import odespy
sys.path.append('../python')
import pyfun as pf

# ------------------------------------------------------------------------------
# Python numerical simulation
# ------------------------------------------------------------------------------
# Parameters
pars = {'s'					: 1.0, # saturation
		'por'				: 0.15, # porosity
		'v_cell'			: 1.0, # m^3_bulk
		'rho_rock'			: 1200.0, # kg/m^3_bulk

		'k_s2o4_disp'		: 0.0, # [/s]
		'k_s2o4_o2'			: 1.e-3, # [/s]
		'k_s2o4_fe3'		: 0.0, # [/s]
		'k_fe2_o2_fast'		: 0.0, # [/s]
		'k_fe2_o2_slow'		: 0.0, # [/s]
		'k_fe2_cr6_fast'	: 0.0, # [/s]
		'k_fe2_cr6_slow'	: 0.0, # [/s]

		'ssa_feoh3'			: 175.0, # [m^2/g]
		'alpha'				: 0.6, # [/s]
		'mv_feoh3'			: 3.436e-5, # [m^3_mnrl/m^3_bulk]
		'mw_feoh3'			: 106.87 # [g/mole]
}

# Solver options
sopt = {'T' :	10.0 * 24 * 60 * 60, # end of simulation [s from d]
		'N' :	1000, # no of timesteps
}

# Initial conditions
init = {'H+'		: 1.e-11, # [M]
		'O2(aq)'	: 1.e-3, # [M]
		'CrO4--'	: 1.e-20, # [M]
		'Cr+++'		: 1.e-20, # [M]
		'S2O4--'	: 1.e-2, # [M]
		'S2O3--'	: 1.e-20, # [M]
		'SO3--'		: 1.e-20, # [M]
		'SO4--'		: 1.e-20, # [M]
		'Fe+++'		: 1.e-20, # [M]
		'Fe++'		: 1.e-20, # [M]
		'fast_Fe++' : 1.e-20, # [m^3/m^3_bulk]
		'slow_Fe++' : 1.e-20, # [m^3/m^3_bulk]
		'Fe(OH)3_s' : 1.e-20, # [m^3/m^3_bulk]
}

dithionite_sandbox = pf.make_dithionite_sandbox(pars)
u, t = pf.run_ode(init, pars, sopt, dithionite_sandbox)

L_water = pars['v_cell'] * pars['por'] * pars['s'] * 1.e3 # [L]
results_ode = {}
results_ode['time'] = t/3600/24 # [d from s]
results_ode['H+'] = u[:,0]/L_water
results_ode['O2(aq)'] = u[:,1]/L_water
results_ode['S2O4--'] = u[:,4]/L_water
results_ode['SO3--'] = u[:,6]/L_water
results_ode['SO4--'] = u[:,7]/L_water

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "s2o4_o2"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total H+ [M]', 'Total S2O4-- [M]', 'Total SO3-- [M]', 'Total SO4-- [M]', 'Total O2(aq) [M]']
observation_list = ['obs (1)']
results_pflotran =  pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# First plot
fig = plt.figure(figsize=[7,5])
skipfactor = 50 # skip data in ode results
xlims = [0,10]
fontsize = 9

# First plot
ax = fig.add_subplot(2, 2, 1)
pflo_plotvars = [[variable_list[1]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['S2O4--']
legend_list = ['S2O4--, PFLOTRAN','S2O4--, odespy']

pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [d]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

regression_result = pf.calc_regression(ts = 1.0,tol = 1.0e-6,results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)

# Second plot
ax = fig.add_subplot(2, 2, 2)
pflo_plotvars = [[variable_list[2]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['SO3--']
legend_list = ['SO3--, PFLOTRAN','SO3--, odespy']

pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [d]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

regression_result = pf.calc_regression(ts = 1.0,tol = 1.0e-6,results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)

# Third plot
ax = fig.add_subplot(2, 2, 3)
pflo_plotvars = [[variable_list[3]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['SO4--']
legend_list = ['SO4--, PFLOTRAN','SO4--, odespy']

pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [d]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

regression_result = pf.calc_regression(ts = 1.0,tol = 1.0e-6,results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)

# Fourth plot
ax = fig.add_subplot(2, 2, 4)
pflo_plotvars = [[variable_list[4]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['O2(aq)']
legend_list = ['O2(aq), PFLOTRAN','O2(aq), odespy']
# 
pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [d]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

regression_result = pf.calc_regression(ts = 1.0,tol = 1.0e-6,results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)

plt.suptitle("S2O4-- oxidized by O2(aq)")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(simbasename + '.png')
