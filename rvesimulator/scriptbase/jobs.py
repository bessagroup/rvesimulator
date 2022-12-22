# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *


class AbaqusJobs(object):   
    """
    Define the jobs of abaqus simulation 
    """
    def __init__(self, model_name, job_name, subtoutine_path=''):
        self.model_name = model_name
        self.job_name = job_name
        self.subroutine_path = subtoutine_path

    def sequential_job(self):
        """
        this function is used to create a sequential job in abaqus 
        """     
        mdb.Job(name=self.job_name, model=self.model_name, description='', type=ANALYSIS, atTime=None,
                waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE,
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE,
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint= OFF,
                contactPrint=OFF, historyPrint=OFF, userSubroutine=self.subroutine_path,
                scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
        mdb.jobs[self.job_name].submit(consistencyChecking=OFF)
        mdb.jobs[self.job_name].waitForCompletion()  