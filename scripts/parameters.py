# Script to set up parameters for SYMPHONY

# SYMPHONY version
versions = ['rvf']

# Output parent path
outputDir = '/Users/feb223/projects/coin/RVF/tests'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
instanceDirs = {
     'KP' : '/Users/feb223/projects/coin/RVF/Data_rvf/KP',
     'SPP' : '/Users/feb223/projects/coin/RVF/Data_rvf/SPP',
     'MILP' : '/Users/feb223/projects/coin/RVF/Data_rvf/MILP',
     'IP' : '/Users/feb223/projects/coin/RVF/Data_rvf/IP'
}

# Set up senarios
# SYMPHONY additional parameters to be set
symParams = {}

symParams['RVFconstruct'] = {
}