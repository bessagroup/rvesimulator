import os
import sys
from collections import OrderedDict

import f3dasm
from matplotlib import pyplot as plt

from rvesimulator.cases.hollow_plate_rve import NaiveHollowPlate

# create the doe
N = 2  # number of samples

# define the doe
radius = f3dasm.ContinuousParameter(
    name="radius", lower_bound=0.1, upper_bound=0.4
)
size = f3dasm.ContinuousParameter(
    name="size", lower_bound=1.0, upper_bound=1.5
)
# define the output
stress = f3dasm.ContinuousParameter(name="stress")
strain = f3dasm.ContinuousParameter(name="strain")

design = f3dasm.DesignSpace(
    input_space=[radius, size], output_space=[stress, strain]
)

sampler = f3dasm.sampling.LatinHypercube(design=design, seed=1)
data = sampler.get_samples(numsamples=N)
print(data.data)

#
simulator = NaiveHollowPlate()
simulator.update_sim_info(strain=[0.05, 0.05, 0.05], print_info=True)


# # calculate the initial responses of simulation
# for ii in range(len(samples_dict)):
#     results = simulator.run_simulation(
#         sample=samples_dict[ii], third_folder_index=ii
#     )
#     # update DoE information
#     for jj in range(len(list(responses.keys()))):
#         responses.at[ii, list(responses.keys())[jj]] = results[
#             list(responses.keys())[jj]
#         ]
results = simulator.run_f3dasm(data=data)
print(results.data)
