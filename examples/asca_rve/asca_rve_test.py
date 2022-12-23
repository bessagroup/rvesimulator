import sys
from collections import OrderedDict

import pandas as pd
from matplotlib import pyplot as plt

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)

# local libraries
import rvesimulator
from examples.asca_rve.simulator_caller import SimulatorCaller
from rvesimulator.design_of_experiment.samplers import FixNumberSampler

# create the doe
# define the design space
doe_variables = OrderedDict({"vol_req": 0.30})

# define number of samples
num_points = 3
# define the information of outputs
name_outputs = ["vol_frac", "PK2", "Green_strain"]
doe_sampler = FixNumberSampler()
doe_sampler.sampling(
    num_samples=num_points,
    design_space=doe_variables,
    out_names=name_outputs,
    seed=123456,
)

print(doe_sampler.data)
#
simulation_wrapper = SimulatorCaller()
simulation_wrapper.update_sim_info(print_info=True)
data = simulation_wrapper.run_simulation(data=doe_sampler.data)
print(data)


samples = data["samples"]
respones = data["responses"]
exp_data = pd.read_csv("Original data.csv", header=None)
pparam = dict(xlabel="$\epsilon_{xx}$", ylabel="$\sigma_{xx}$ (MPa)")
with plt.style.context(["science", "ieee"]):
    fig, ax = plt.subplots()
    ax.plot(exp_data.loc[:, 0], exp_data.loc[:, 1], label="Reference results")
    ax.plot(
        respones["Green_strain"][0][:, 0, 0],
        respones["PK2"][0][:, 0, 0],
        label="realization 1",
    )
    ax.plot(
        respones["Green_strain"][1][:, 0, 0],
        respones["PK2"][1][:, 0, 0],
        label="realization 2",
    )
    ax.plot(
        respones["Green_strain"][2][:, 0, 0],
        respones["PK2"][2][:, 0, 0],
        label="realization 2",
    )
    ax.legend()
    ax.set(**pparam)
    fig.savefig("fig_sigma_xx.pdf")
    fig.savefig("fig_sigma_xx.png", dpi=300)
