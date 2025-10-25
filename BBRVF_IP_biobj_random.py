#!/usr/bin/env python
# coding: utf-8

# In[50]:


'''
Created on Dec 12, 2023

@author: Samira Fallah (saf418@lehigh.edu)
'''


# In[51]:


import time
import argparse
from gurobipy import *
import math, sys, itertools, os

# Absolute path to the libSym, extention change based on the OS:
# - Windows: .dll
# - Linux: .so
# - OSX: .dylib


# In[52]:


# This is the path pointing to the dynamic library
# lib_path = r'/mnt/c/rvf_symphony/build-SYMPHONY-rvf/lib/libSym.so'
# server
lib_path = r'/Users/feb223/projects/coin/RVF/build-SYM-opt/lib/libSym.dylib'

# Import the python package
from symphony import *

# Load the library by calling this static method
Symphony.dlopen(lib_path)


# server
# sys.path.insert(0, r'/home/saf418/ValueFunctionCode/RVFCodes/SYMPHONYBBRVF/Test_py_IP_biobj')

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--instance", help = "specify the instance")
parser.add_argument("-p", "--policy", help = "policy for tree selection", default='0')
parser.add_argument("-n", "--samples", help = "number of samples for zeta", default=20)
parser.add_argument("-t", "--timelimit", help = "time limit in seconds", default=1e50)
flags = parser.parse_args()
instance = flags.instance 
tree_policy = flags.policy
num_samples = int(flags.samples)
timeLimit = int(flags.timelimit)

sys.path.insert(0, os.path.dirname(instance))

instanceName = (instance.split("/")[-1]).split('.mps')[0]
for key, val in vars(__import__(instanceName)).items():
    if key.startswith('__') and key.endswith('__'):
        continue
    vars()[key] = val

model_path = instance

debug_print = False
ipForU = True


def changeValue(value):
    if str(value) == 'None':
        return 0.0
    return value

# Generate a feasible solution to start with
def generateInitPoint():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.INTEGER, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS), GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
   
    temp = {}
    for k, v in intVarsInit.items():
        temp[k] = changeValue(v.X) 
    
    return temp



# Generate U - Upper bound on the whole VF
def generateU():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    if ipForU:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
        slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.INTEGER, lb = 0, name = "slack variable")
    else:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.CONTINUOUS, lb = 0, ub = UB_I, name = "integer variable")
        slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
        
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS), GRB.MAXIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
    
    return m.objVal




# Generate L - Lower bound on the whole VF
def generateL():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    if ipForU:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
        slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.INTEGER, lb = 0, name = "slack variable")
    else:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.CONTINUOUS, lb = 0, ub = UB_I, name = "integer variable")
        slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
        
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS), GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
    
    return m.objVal




def generateLBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS), GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()

    return m.objVal



def generateUBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS), GRB.MAXIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
    
    return m.objVal




# Convert the feasible solution to an efficient solution
def convertWeakToStrongNDP(_intVarsInit, _print=False):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsStrong = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    slackVarsStrong = m.addVars(list(SLACKVARS), vtype = GRB.INTEGER, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsStrong[i] for i in INTVARS) 
                   + sum(MAT[(i, j)] * intVarsStrong[j] for j in INTVARS for i in CONSVARRHS), GRB.MINIMIZE)
    
    RHSFirstObj = 0
    for i in INTVARS:
        RHSFirstObj += OBJ[i] * _intVarsInit[i]
      
    m.addConstr(sum(OBJ[i] * intVarsStrong[i] for i in INTVARS) <= RHSFirstObj)
    
    for j in CONSVARRHS:
        m.addConstr(sum(MAT[(j, i)] * intVarsStrong[i] for i in INTVARS) <= sum(MAT[(j, i)] * _intVarsInit[i] for i in INTVARS))
        
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsStrong[i] for i in INTVARS)  
                    + sum(MATFixed[(j, i)] * slackVarsStrong[i] for i in SLACKVARS)  == RHS[j])
    
    m.optimize()
    
    temp = {}
    for k, v in intVarsStrong.items():
        temp[k] = changeValue(v.X)
    
    return temp



# enhancement idea 2 of the linear version
def RVFSubproblem(timelimit=None):

    m = Model()
    m.setParam("LogToConsole", 0);
    if timelimit:
        m.setParam("TimeLimit", timelimit);
    # This parameter is redundant; you can remove it. 
    # Since the objective is an arbitrary value, the first feasible solution would be the optimal one.
    m.params.SolutionLimit = 1
    
    thetaVar = m.addVar(vtype = GRB.CONTINUOUS, name = "theta")
    intVars = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "int_var")
    slackVars = m.addVars(list(SLACKVARS), vtype = GRB.INTEGER, lb = 0, name = "slack_var")
    alphaVars = m.addVars(len(CONSVARRHS), len(intPartList), vtype = GRB.BINARY, name = "alpha_var")
    betaVars = m.addVars(len(intPartList), vtype = GRB.BINARY, name = "beta_var")
    
    m.setObjective(1, GRB.MAXIMIZE)
    
    # theta constraint -- const 1
    for k in range(len(intPartList)):
        m.addConstr(thetaVar <= - sum(OBJ[i]*intVars[i] for i in INTVARS) +   
            (1 - betaVars[k])*sum(OBJ[j]*intPartList[k][j] for j in INTVARS) + betaVars[k]*(U+1))
    
    # const 2 
    for k in range(len(intPartList)):
        for i in CONSVARRHS:
            m.addConstr(sum(MAT[(i, j)]*(intVars[j] - intPartList[k][j]) for j in INTVARS) + eps 
                                                    <= M[i]*(1 - alphaVars[(i, k)]))
    
    # const 3 
    for k in range(len(intPartList)):
        for i in CONSVARRHS:
            m.addConstr(sum(MAT[(i, j)]*(intPartList[k][j] - intVars[j]) for j in INTVARS) 
                            <= M[i]*alphaVars[(i, k)])
    
    # const 4
    for k in range(len(intPartList)):
        for i in CONSVARRHS:
            m.addConstr(betaVars[k] >= alphaVars[(i, k)])
    
    # const 5
    for k in range(len(intPartList)):
        m.addConstr(betaVars[k] <= sum(alphaVars[(i, k)] for i in CONSVARRHS))
    
    # const 6
    for k in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(k, j)]*intVars[j] for j in INTVARS) + 
                    sum(MATFixed[(k, j)]*slackVars[j] for j in SLACKVARS) == RHS[k])

    m.addConstr(thetaVar >= 0.90)
    
    m.optimize()

    status = m.getAttr('Status')
    
    if status in [3, 4]:
        return -1
    
    elif status == 9:
        return -2

    temp_ndp = () 

    int_part = convertWeakToStrongNDP(dict((k, round(changeValue(v.X))) for k,v in intVars.items()), _print=False)
    
    temp_ndp = temp_ndp + (round(sum(OBJ[j]*int_part[j] for j in INTVARS)),)
        
    for k in CONSVARRHS:
        temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*int_part[l] for l in INTVARS)),) 
    
    # print('New NDP (RVF): ', temp_ndp)
    EF.append(temp_ndp)    
    intPartList.append(int_part)
    
    return None




fixed_RHS_list = list(RHS.values())

eps = 0.5

M = {}

for i in CONSVARRHS:
    M[i] = sum(abs(MAT[(i, j)]) for j in range(numIntVars)) + 1

RHS_list = []
RHS_list.append(0) 

for j in range(len(fixed_RHS_list)):
    RHS_list.append(fixed_RHS_list[j])



def calcRVFUBAtZeta(sampleZeta):
    filteredEF = [item for item in EF if item[1] <= sampleZeta]
    # print('filteredEF', filteredEF)
    # print('EF', EF)
    if not filteredEF:
        return U
    else:
        return min([item[0] for item in filteredEF])




def findGapOfZeta(sampleZeta):
    
    UB_bound = calcRVFUBAtZeta(sampleZeta[0])

    st = time.time()

    LB_bound = forest.evaluate_forest_dual_function(sampleZeta)

    et = time.time()

    elapsed = et - st
    global totTimeDfEval

    totTimeDfEval += elapsed

    if LB_bound >= 1e+20:
        return UB_bound, None
            
    return UB_bound, LB_bound


# In[77]:


def generate_sequence(lst):
    result = []
    result.extend([lst[0], lst[-1]])
    result.append((lst[0] + lst[-1]) / 2)
    result = sorted(result)

    for i in range(2):
        temp_lst = []
        for j in range(len(result)-1):
            mid_value = (result[j] + result[j+1]) / 2
            temp_lst.append(mid_value)
            
        result.extend(temp_lst)
        
        result = sorted(result)
        
    return result


minZeta = generateLBZeta(0)
maxZeta = generateUBZeta(0)

print('minZeta:', minZeta)
print('maxZeta:', maxZeta)

my_list = []
my_list.append(minZeta)
my_list.append(maxZeta)
tempListZeta = []


tempListZeta = generate_samples_hypercube(num_samples, [minZeta], [maxZeta], 1, method='lhs')


sym = None
forest = SymphonyForest()
forest.model = model_path

if tree_policy == '0':
    # Always warm start
    forest.params["max_opt_gap"] = INF
    forest.params["max_tree_size"] = INF
elif tree_policy == '1':
    # Never warm start
    forest.params["max_opt_gap"] = 0.00001
    forest.params["max_tree_size"] = INF
elif tree_policy == '2':
    # Hybrid
    forest.params["max_opt_gap"] = 10
    forest.params["max_tree_size"] = 1000
else:
    print('Unknown tree selection policy!')
    exit(1)

totalTime = 0
totTimeDfEval = 0
totTimeSelectTree = 0
remainingTime = 0
intPartList = []
EF = []
mainIter = 0
RVFIter = 0 

U = generateU()
L = generateL()

hitTimeLimit = False
st = time.time()

sym = forest.add_tree()

et = time.time()
timeElapsed = et - st 
remainingTime = timeLimit - timeElapsed
sym.set_param("time_limit", remainingTime)

if (sym.solve() == FUNCTION_TERMINATED_ABNORMALLY):
    print("Something went wrong!")
    # exit(1)  
    
sym.build_dual_function()

opt_sol = sym.get_col_solution()
if opt_sol:
    int_part = convertWeakToStrongNDP(opt_sol, _print=False)
    temp_ndp = ()
    temp_ndp = temp_ndp + (round(sum(OBJ[j]*int_part[j] for j in INTVARS)),)

    for k in CONSVARRHS: 
        temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*int_part[l] for l in INTVARS)),)

    if temp_ndp not in EF:
        # print('New NDP (SYMPHONY): ', temp_ndp)
        EF.append(temp_ndp)
        intPartList.append(int_part)



feas_status = 0
EF_len_main = 0
EF_len_after_main = 0
numLBGRVF = 0
LBGRVFLst = []

for sample_zeta in tempListZeta:
    
    mainIter += 1
    print('------------------------')
    # print('iteration', iter)
    print('sample_zeta', sample_zeta)
    # Convert to a python list 
    sample_zeta = sample_zeta.astype(float).tolist()
    UB, LB = findGapOfZeta(sample_zeta)
    
    if LB == None:
        continue
        
    gap_zeta = UB - LB
            
    if (gap_zeta <= 0.1):
        RVFIter += 1
        et = time.time()
        timeElapsed = et - st 
        remainingTime = timeLimit - timeElapsed
        feas_status = RVFSubproblem(timelimit=remainingTime)
        # print('RVFSubproblem Called')
        if feas_status == -1:
            print('Reached to optimality!')
            break
        elif feas_status == -2:
            print('Time limit surpassed!')
            break
        else:
            continue
       
    select_tree_start = time.time()
    sym = forest.pick_tree(sample_zeta, UB)
    select_tree_end = time.time()
    elapsed_select_tree = select_tree_end - select_tree_start
    totTimeSelectTree += elapsed_select_tree

    for i in CONSVARRHS:
        sym.set_row_upper(i, sample_zeta[i])

    et = time.time()
    timeElapsed = et - st 
    remainingTime = timeLimit - timeElapsed
    sym.set_param("time_limit", remainingTime)
        
    if (sym.warm_solve() == FUNCTION_TERMINATED_ABNORMALLY):
        print("Something went wrong!")
        # exit(1)

    sym.build_dual_function()
        
    opt_sol = sym.get_col_solution()
    
    if opt_sol:
        int_part = convertWeakToStrongNDP(opt_sol, _print=False)
        temp_ndp = ()
        temp_ndp = temp_ndp + (round(sum(OBJ[j]*int_part[j] for j in INTVARS)),)
    
        for k in CONSVARRHS: 
            temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*int_part[l] for l in INTVARS)),)

        if temp_ndp not in EF:
            # print('New NDP (SYMPHONY): ', temp_ndp)
            EF.append(temp_ndp)
            intPartList.append(int_part)

    EF_len_main = len(EF) 

    # Check time limit
    et = time.time()
    timeElapsed = et - st 
    if timeElapsed > timeLimit:
        hitTimeLimit = True
        print('Time limit surpassed!')
        break
    
print('---------------------------------------------')
print('Main loop ended')
print('---------------------------------------------')

# print('RVF is called from now on')
while feas_status != -1 and not hitTimeLimit:
    # totalIter += 1
    RVFIter += 1
    # print('iteration', iter)
    et = time.time()
    timeElapsed = et - st 
    remainingTime = timeLimit - timeElapsed
    feas_status = RVFSubproblem(timelimit=remainingTime)
    # print('RVFSubproblem Called')
    # print('------------------------')
    if feas_status == -1:
        print('Reached to optimality!')
        break
    elif feas_status == -2:
        print('Time limit surpassed!')
        break

    # Check time limit
    et = time.time()
    timeElapsed = et - st 
    if timeElapsed > timeLimit:
        hitTimeLimit = True
        print('Time limit surpassed!')
        break
    
        
EF_len_after_main = len(EF) -  EF_len_main


et = time.time()
totalTime += et - st 
totalTimeWOLP = 0




print('---------------------------------------------')
print('Problem structure')
print('---------------------------------------------')
print('       Total number of vars: ', numVars)
print('                    Integer: ', numIntVars)
print('                 Continuous: ', 0)
print('Total number of constraints: ', numConsFixed)
print(' Total number of objectives: ', numConsVar + 1)

print('---------------------------------------------')
print('Statistics')
print('---------------------------------------------')


print('Total Running Time (sec): ', totalTime)
print('Total Running Time w/o LPs (sec): ', totalTimeWOLP)
print('EF: ' , EF)
print('Number of total NDPs: ', len(EF))
print('Number of main algo iterations: ', mainIter)
print('Number of iterations that RVF called: ', RVFIter)
print('Number of total iterations: ', mainIter + RVFIter)
print('Len EF inside the main: ', EF_len_main)
print('Len EF after the main: ', EF_len_after_main)
print('Number of cases that LB is greater than RVF: ', numLBGRVF)
print('List of zetas that LB is greater than RVF: ', LBGRVFLst)
print('---------------------------------------------')
print('Forest Statistics')
print('---------------------------------------------')
print('Total number of trees: ', len(forest.trees))
print('Total time evaluating dual functions (sec): ', totTimeDfEval)
print('Total time selecting trees (sec): ', totTimeSelectTree)

# with open("Results_IP_biobj_details_{}.txt".format(instance.split('.')[0]), "a") as _filed:
#     _filed.write('Total Running Time (sec): ' + str(totalTime) + '\n') 
#     _filed.write('Total Running Time w/o LPs (sec): ' + str(totalTimeWOLP) + '\n') 
#     _filed.write('EF: ' + str(EF) + '\n') 
#     _filed.write('Number of total NDPs: ' + str(len(EF)) + '\n') 
#     _filed.write('Number of main algo iterations: ' + str(mainIter) + '\n') 
#     _filed.write('Number of iterations that RVF called: ' + str(RVFIter) + '\n') 
#     _filed.write('Number of total iterations: ' + str(mainIter + RVFIter) + '\n') 
#     _filed.write('Len EF inside the main: ' + str(EF_len_main) + '\n') 
#     _filed.write('Len EF after the main: ' + str(EF_len_after_main) + '\n') 
#     _filed.write('Number of cases that LB is greater than RVF: ' + str(numLBGRVF) + '\n') 
#     _filed.write('List of zetas that LB is greater than RVF: ' + str(LBGRVFLst)) 



del sym
