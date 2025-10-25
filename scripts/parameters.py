# Script to set up parameters for SYMPHONY

# SYMPHONY version
versions = ['forest']

# Output parent path
outputDir = '/Users/feb223/projects/coin/RVF/tests'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
instanceDirs = {
     # 'KP' : '/Users/feb223/projects/coin/RVF/Data_rvf/KP',
     # 'SPP' : '/Users/feb223/projects/coin/RVF/Data_rvf/SPP',
     # 'MILP' : '/Users/feb223/projects/coin/RVF/Data_rvf/MILP',
     # 'IP' : '/Users/feb223/projects/coin/RVF/Data_rvf/IP'
     'Sample' : '/Users/feb223/projects/coin/RVF/pySYMPHONY/Data_rvf/Sample'
}

# Set up senarios
commonParams = {
     'timelimit' : '1800'
}
# SYMPHONY additional parameters to be set
symParams = {
    
}

symParams['forest_warmstart'] = {
    'policy' : '0'
}

symParams['forest_coldstart'] = {
    'policy' : '1'
}

symParams['forest_hybrid'] = {
    'policy' : '2'
}