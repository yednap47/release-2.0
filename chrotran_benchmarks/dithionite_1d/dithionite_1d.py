import sys
sys.path.append('../python')
from matplotlib import pyplot as plt
import pyfun as pf
import numpy as np
import itertools as it

simbasename = "dithionite_1d"
observation_filename = [simbasename + '-mas.dat']
variable_list = ["east CrO4-- [mol/d]"]
observation_list = ['']
results_pflotran = pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

gold_name = [simbasename + '-mas.dat.gold']
gold_standard = pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=gold_name)

fig = plt.figure(figsize=[5,4])
majorFormatter = plt.matplotlib.ticker.FormatStrFormatter('%0.2e')
mycmap=plt.cm.jet(np.linspace(0,1,5))

plt.plot(results_pflotran['time'], map(lambda x: -x, results_pflotran['east CrO4-- [mol/d] ']),ls = '-',c = mycmap[0], label = "new simulation")
plt.scatter(gold_standard['time'][1:-1:10], map(lambda x: -x, gold_standard['east CrO4-- [mol/d] '][1:-1:10]),label = "gold standard")
plt.xlim(0,365)
plt.xlabel("Time [d]")
plt.ylabel(variable_list[0])
ax = plt.gca()
ax.yaxis.set_major_formatter(majorFormatter)
ax.legend(ncol=1, fancybox=True, shadow=False, prop={'size': '10.0'}, loc='best')
plt.tight_layout()

plt.savefig('dithionite_' + simbasename + '.png',dpi=100)

# For regression test, use gold standard file instead of python ODE
pflo_plotvars = [[variable_list[0]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars  = [variable_list[0] + " "]
regression_result = pf.calc_regression(ts = 20.0,tol = 1.0e-9,results_ode=gold_standard, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, sim=simbasename)
