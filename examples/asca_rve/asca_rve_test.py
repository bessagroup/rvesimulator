from collections import OrderedDict

import pandas as pd
from femai.DesignOfExperiment.Samplers import FixNumberSampler
from matplotlib import pyplot as plt

## local libraries
from examples.asca_rve.simulator_caller import SimulatorCaller

## create the doe
# define the design space
doe_variables = OrderedDict({'Volume_req': 0.30})

# define number of samples
num_points = 3
# define the information of outputs
name_outputs = ['Volume_frac', 'PK2', 'Green_strain']
doe_sampler = FixNumberSampler()
samples = doe_sampler.Sampling(num_samples=num_points,
                               design_space=doe_variables,
                               out_names=name_outputs)
print(samples)
#
simulation_wrapper = SimulatorCaller()
simulation_wrapper.UpdateSimulationInformation(print_info=True)
samples = simulation_wrapper.RunSimulation(DoE=samples)
print(samples)

exp_data = pd.read_csv('Original data.csv', header=None)
pparam = dict(xlabel='$\epsilon_{xx}$', ylabel='$\sigma_{xx}$ (MPa)')
with plt.style.context(['science', 'ieee']):
    fig, ax = plt.subplots()
    ax.plot(exp_data.loc[:, 0], exp_data.loc[:, 1], label='Reference results')
    ax.plot(samples['Green_strain'][0][:, 0, 0],  samples['PK2'][0][:, 0, 0],  label='point_0')
    ax.plot(samples['Green_strain'][1][:, 0, 0],  samples['PK2'][1][:, 0, 0],  label='point_1')
    ax.plot(samples['Green_strain'][2][:, 0, 0],  samples['PK2'][2][:, 0, 0],  label='point_2')
    ax.legend()
    ax.set(**pparam)
    fig.savefig('fig1.pdf')
    fig.savefig('fig1.png', dpi=600)



