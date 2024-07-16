# abaqus
from abaqus import *
from abaqusConstants import *
from caeModules import *


class Jobs:
    def create_sequential_job(self, subroutine_path):
        """create sequential jobs

        Parameters
        ----------
        subroutine_path : str
            path for subroutine
        """

        mdb.Job(
            name=self.job_name,
            model=self.model_name,
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine=subroutine_path,
            scratch="",
            resultsFormat=ODB,
            multiprocessingMode=DEFAULT,
            numCpus=self.num_cpu,
            numGPUs=0,
            numDomains=self.num_cpu)
        mdb.jobs[self.job_name].writeInput(consistencyChecking=OFF)

    def submit_job(self):
        mdb.jobs[self.job_name].submit(consistencyChecking=OFF)
        mdb.jobs[self.job_name].waitForCompletion()
