abaqus2py
=========

This module is most important part of this repo, it used to run abaqus python script with out abaqus GUI.
Ideally, this module can handle any abaqus simulation as long as you have a well defined abaqus script. 
However, in order to make your script work with this module, you need to modify your script a little bit.

script requirements
-------------------

1. if you prefer to use **functional** code, you have to modify your script into the following format:

.. code-block:: python
    :name: abaqus_sim_function.py

    
    def simulation(sim_info):
        # note: name of the function is optinal, you can name it whatever you want
        # sim_info is a dictionary contains all variables you need to run your simulation
        # for example:
        # sim_info = {"youngs_modulus": 1000, "poisson_ratio": 0.3}
        youngs_modulus = sim_info["youngs_modulus"]
        poisson_ratio = sim_info["poisson_ratio"]
        # write your whole abaqus script below:
        ...
    
    def post_process(job_name):
        # note: name of the function is optinal, you can name it whatever you want
        # job_name is the name of the job you just run

        # write your post process code here
        ...

        # finally, return the result you want to save into pickle file
        results = {"stress": stress, "strain": strain}
        with open("results.p", "w") as f:
            pickle.dump(results, f)

2. if you prefer to use **objected-oriented** code, you have to modify your script into the following format:

.. code-block:: python 
    :name: abaqus_sim_main.py


    class Simulation(object):
        def __init__(self, sim_info):
            # sim_info is a dictionary contains all variables you need to run your simulation
            # for example:
            # sim_info = {"youngs_modulus": 1000,
            #              "poisson_ratio": 0.3,
            #              "platform": "ubuntu",
            #              "job_name": "job"}
            self.youngs_modulus = sim_info["youngs_modulus"]
            self.poisson_ratio = sim_info["poisson_ratio"]
            self.platform = sim_info["platform"]
            self.job_name = sim_info["job_name"]
            # write write other variables you need here
            ...
            # do simulation
            self.simulation()
            # do raw post process right after simulation
            if self.platform == "cluster":
                PostProcess(self.job_name)


        def simulation(self):
            # write your whole abaqus script below:
            ...

    class PostProcess(object):
        def post_process(self, job_name):
            # job_name is the name of the job you just run
            # write your post process code here
            ...

            # finally, return the result you want to save into pickle file
            results = {"stress": stress, "strain": strain}
            with open("results.p", "w") as f:
                pickle.dump(results, f)

basic usage
-----------
Once you have your script ready as the following format, then you can put your script to a folder you want. Here I want 
to show a simple example about the folder structure of your own script and :any:`abaqus2py` module.

.. literalinclude:: abaqus2py.txt
   :language: text
   :lines: 1-13

The above case, the folder **abaqus_script** is parallel to the repo folder **rvesimulator**. 
Then you can run your simulation by:

.. code-block:: python

    from rvesimulator.abaqus2py.abaqus_simulator import AbaqusSimulator
    folder_info = {
    "main_work_directory": os.path.join(os.getcwd(), "Data"), # you can change where to do your simulation here
    "script_path": ".../abaqus_script", # absolute path to your script folder
    "current_work_directory": "data", # random name
    "sim_path": "case", # important part, this is the name of your script
    "sim_script": "Simulation", # important
    "post_path": "case",# important
    "post_script": "PostProcess",# important
    }
    sim_info = {
                job_name: "job", # important
                platform: "ubuntu", # important
                }
    # initiate the simulator
    simulator = AbaqusSimulator(sim_info=sim_info, folder_info=folder_info)
    # run the simulation
    simulator.run()
    # read back the results
    results = simulator.read_back_results()

.. note::
    where ``folder_info`` is a dictionary indicates where to find your abaqus script and where to run and save your simulation.
    In the way, it has a bit strict requirement on the folder structure of your script. For the ``sim_info`` dictionary, it encode the design 
    of experiment to abaqus simulation. Also, the **job_name** and **platform** are important keys in the dictionary, you have to include them in your 
    ``sim_info`` dictionary.



    


