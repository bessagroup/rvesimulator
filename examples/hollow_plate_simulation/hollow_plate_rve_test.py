import sys
from collections import OrderedDict

from matplotlib import pyplot as plt

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)

from rvesimulator.design_of_experiment.samplers import RandomSampler
# from simulator_caller import SimulatorCaller
from rvesimulator.simulators.hollow_plate_rve import HollowPlateRVE

# create the doe
# define the design space
doe_variables = OrderedDict({"radius": [0.10, 0.40], "size": [1.0, 1.5]})
# define number of samples
num_points = 3
# define the information of outputs
name_outputs = ["PK2", "Green_strain"]
doe_sampler = RandomSampler()
doe_sampler.sampling(
    num_samples=num_points,
    design_space=doe_variables,
    out_names=name_outputs,
    seed=123456,
)

print(doe_sampler.data)
#
simulation_wrapper = HollowPlateRVE()
simulation_wrapper.update_sim_info(loads=[0.05, 0.05, 0.05])

data = simulation_wrapper.run_simulation(data=doe_sampler.data)

print(data)
samples = data["samples"]
respones = data["responses"]
# define parameter of plot figures
pparam = dict(xlabel="$\epsilon_{xx}$", ylabel="$\sigma_{xx}$ (MPa)")
with plt.style.context(["science", "ieee"]):
    fig, ax = plt.subplots()
    ax.plot(
        respones["Green_strain"][0][:, 0, 0],
        respones["PK2"][0][:, 0, 0],
        label=f'radius= {samples["radius"][0]:3f}',
    )
    ax.plot(
        respones["Green_strain"][1][:, 0, 0],
        respones["PK2"][1][:, 0, 0],
        label=f'radius= {samples["radius"][1]:3f}',
    )
    ax.plot(
        respones["Green_strain"][2][:, 0, 0],
        respones["PK2"][2][:, 0, 0],
        label=f'radius= {samples["radius"][2]:3f}',
    )
    ax.legend()
    ax.set(**pparam)
    fig.savefig("fig_sigma_xx.pdf")
    fig.savefig("fig_sigma_xx.png", dpi=300)

with plt.style.context(["science", "ieee"]):
    fig, ax = plt.subplots()
    ax.plot(
        respones["Green_strain"][0][:, 1, 1],
        respones["PK2"][0][:, 1, 1],
        label=f'radius= {samples["radius"][0]:3f}',
    )
    ax.plot(
        respones["Green_strain"][1][:, 1, 1],
        respones["PK2"][1][:, 1, 1],
        label=f'radius= {samples["radius"][1]:3f}',
    )
    ax.plot(
        respones["Green_strain"][2][:, 1, 1],
        respones["PK2"][2][:, 1, 1],
        label=f'radius= {samples["radius"][2]:3f}',
    )
    ax.legend()
    ax.set(**pparam)
    fig.savefig("fig_sigma_yy.pdf")
    fig.savefig("fig_sigma_yy.png", dpi=300)


with plt.style.context(["science", "ieee"]):
    fig, ax = plt.subplots()
    ax.plot(
        respones["Green_strain"][0][:, 1, 0],
        respones["PK2"][0][:, 1, 0],
        label=f'radius= {samples["radius"][0]:2f}',
    )
    ax.plot(
        respones["Green_strain"][1][:, 1, 0],
        respones["PK2"][1][:, 1, 0],
        label=f'radius= {samples["radius"][1]:2f}',
    )
    ax.plot(
        respones["Green_strain"][2][:, 1, 0],
        respones["PK2"][2][:, 1, 0],
        label=f'radius= {samples["radius"][2]:2f}',
    )
    ax.legend()
    ax.set(**pparam)
    fig.savefig("fig_sigma_xy.pdf")
    fig.savefig("fig_sigma_xy.png", dpi=300)
