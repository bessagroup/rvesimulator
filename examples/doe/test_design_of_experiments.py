# third-party packages
import sys
from collections import OrderedDict

# path of local project
folder_path = "/home/jiaxiangyi/Documents/rvesimulator"
sys.path.insert(0, folder_path)
# local packages
from rvesimulator.design_of_experiment.samplers import SobolSequence

# define the design space
doe_variables = OrderedDict({"x1": [0.0, 1.0], "x2": [1.0, 2.0]})
# define number of samples
num_points = 10
# define the information of outputs
name_outputs = ["y"]
doe_sampler = SobolSequence()
doe_sampler.sampling(
    num_samples=num_points,
    design_space=doe_variables,
    out_names=name_outputs,
    seed=123456,
)
print(doe_sampler.data)

doe_sampler.plot_samples(fig_name="sobol", save_fig=True)
