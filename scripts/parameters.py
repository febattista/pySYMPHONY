# Script to set up parameters for SYMPHONY

# SYMPHONY version
versions = ['rvf']

# Output parent path
outputDir = '/home/feb223/rvf/tests'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
instanceDirs = {
     'KP' : '/home/feb223/rvf/Data_rvf/KP_MPS',
     'SPP' : '/home/feb223/rvf/Data_rvf/SPP_MPS',
}

# Set up senarios
# SYMPHONY additional parameters to be set
symParams = {}

symParams['dualFunc'] = {
}