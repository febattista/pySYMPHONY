import sys, os, collections
import shutil, subprocess

sys.path.append("..")

import numpy as np
from parameters import versions, outputDir, instanceDirs, symParams

def runExperiments(exe, instanceDirs, outputDir, versions, symParams):
    # set up output directories
    # use hierarchy:  outDir/version/param_scenario_name/testset_name/
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    for v in versions:
        currpath = os.path.join(outputDir, v)
        if not os.path.exists(currpath):
            os.mkdir(currpath)
        for scenario in symParams:
            print(scenario)
            currpath = os.path.join(outputDir, v, scenario)
            if not os.path.exists(currpath):
                os.mkdir(currpath)
            for testset in instanceDirs:
                currsubpath1 = os.path.join(currpath, testset)
                if not os.path.exists(currsubpath1):
                    os.mkdir(currsubpath1)

    # Run experiments
    for v in versions:
        for testset in instanceDirs:
            for scenario in symParams:

                paramcmd = " ".join([f"-{key} {value}" for key, value in symParams[scenario].items()])
                print(paramcmd)
                outsubpath = os.path.join(outputDir, v, scenario, testset)
                print(outsubpath)
                with os.scandir(instanceDirs[testset]) as inst_it: 
                    for instance in inst_it:
                        # print(instance.name)
                        if instance.name.lower().endswith('.mps'):
                            outname = instance.name[:-4]+'.out'
                            if any(x in instance.name.lower() for x in ('spp', 'mis')):
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_IP_biobj.py"
                            elif 'kp' in instance.name.lower():
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_IP_multiobj.py"
                            elif 'ip_obj_2' in instance.name.lower():
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_IP_biobj_random.py"
                            elif 'ip_obj_' in instance.name.lower():
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_IP_multiobj_random.py"
                            elif 'obj_2' in instance.name.lower():
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_MILP_biobj.py"
                            else:
                                exe = "/Users/feb223/projects/coin/RVF/pySYMPHONY/BBRVF_MILP_multiobj.py"
                            command = ['python', exe, "-i", instance.path] + paramcmd.split()
                            outfile = open(os.path.join(outsubpath, outname), 'w')
                            try:
                                print('Running {} ...'.format(instance.name + " " + exe))
                                # Run the child command with subprocess
                                result = subprocess.run(command, stdout=outfile, 
                                                        stderr=subprocess.PIPE, check=True)
                                return_code = result.returncode
                            except subprocess.CalledProcessError as e:
                                # If the child command returns a non-zero exit code
                                print("Child command failed with return code:", e.returncode)
                                print("Error output:", e.stderr.decode())
                                return_code = e.returncode
                            outfile.close()

                            # Check the return code
                            if return_code == 139:  # Segmentation fault (SIGSEGV)
                                print("Child command crashed with a segmentation fault.")
                            elif return_code != 0:  # Any other non-zero return code
                                print("Child command failed with return code:", return_code)
                            else:
                                print("Child command executed successfully.")
                            
                            print('Complete {}'.format(instance.name))

if __name__ == "__main__":

    exe = ''
    runExperiments(exe, instanceDirs, outputDir, versions, symParams)