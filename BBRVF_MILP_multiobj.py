#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
Created on Dec 12, 2023

@author: Samira Fallah (saf418@lehigh.edu)
'''


# In[2]:


import time
import argparse
from gurobipy import *
import math, sys, itertools, os

# Absolute path to the libSym, extention change based on the OS:
# - Windows: .dll
# - Linux: .so
# - OSX: .dylib


# In[3]:


# This is the path pointing to the dynamic library
lib_path = r'/Users/feb223/projects/coin/RVF/build-SYM-opt/lib/libSym.dylib'
# lib_path = r'/Users/feb223/projects/coin/RVF/build-SYMPHONY-rvf/lib/libOsiSym.dylib'
# server
# lib_path = r'/home/saf418/ValueFunctionCode/RVFCodes/SYMPHONYBBRVF/build-SYMPHONY-rvf/lib/libSym.so'

# Import the python package
from symphony import *

# Load the library by calling this static method
Symphony.dlopen(lib_path)


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
mipForU = True

defaultTheta = 0.1
currentTheta = 0.1


# In[11]:


def changeValue(value):
    if str(value) == 'None':
        return 0.0
    return value


# In[12]:


# Generate a feasible solution to start with
def generateInitPoint():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS) + sum(OBJ[i] * contVarsInit[i] for i in CONVARS),
                                                           GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
#     print(m.getAttr('Status'))
    
    temp = {}
    for k, v in intVarsInit.items():
        temp[k] = changeValue(v.X) 
    for k, v in contVarsInit.items():
        temp[k] = changeValue(v.X) 
    
    return temp


# In[13]:


# Generate U - Upper bound on the whole VF
def generateU():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    if mipForU:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    else:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.CONTINUOUS, lb = 0, ub = UB_I, name = "integer variable")

    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS) + sum(OBJ[i] * contVarsInit[i] for i in CONVARS),
                                                           GRB.MAXIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
    
    return math.ceil(m.objVal)


# In[14]:


# Generate L - Upper bound on the whole VF
def generateL():
    m = Model()
    m.setParam("LogToConsole", 0);
    
    if mipForU:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    else:
        intVarsInit = m.addVars(list(INTVARS), vtype = GRB.CONTINUOUS, lb = 0, ub = UB_I, name = "integer variable")

    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsInit[i] for i in INTVARS) + sum(OBJ[i] * contVarsInit[i] for i in CONVARS),
                                                           GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()
    
    return math.ceil(m.objVal)


# In[15]:


def generateLBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS) + sum(MAT[(objInd, i)] * contVarsInit[i] for i in CONVARS)
                                                                           , GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()

    return m.objVal


# In[16]:


def generateUBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS) + sum(MAT[(objInd, i)] * contVarsInit[i] for i in CONVARS)
                                                                           , GRB.MAXIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()

    return m.objVal


# In[17]:


# Convert the feasible solution to an efficient solution
def convertWeakToStrongNDP(_totalVars, _print=False):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsStrong = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    contVarsStrong = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsStrong = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * intVarsStrong[i] for i in INTVARS) + sum(OBJ[i] * contVarsStrong[i] for i in CONVARS) 
                   + sum(MAT[(i, j)] * intVarsStrong[j] for j in INTVARS for i in CONSVARRHS)
                   + sum(MAT[(i, j)] * contVarsStrong[j] for j in CONVARS for i in CONSVARRHS), GRB.MINIMIZE)
    
    RHSFirstObj = 0
    for i in INTVARS:
        RHSFirstObj += OBJ[i] * _totalVars[i]
    for i in CONVARS:
        RHSFirstObj += OBJ[i] * _totalVars[i]
      
    m.addConstr(sum(OBJ[i] * intVarsStrong[i] for i in INTVARS) + sum(OBJ[i] * contVarsStrong[i] for i in CONVARS)
                                                            <= RHSFirstObj)
    
    for j in CONSVARRHS:
        m.addConstr(sum(MAT[(j, i)] * intVarsStrong[i] for i in INTVARS) + 
                    sum(MAT[(j, i)] * contVarsStrong[i] for i in CONVARS) <= 
                    sum(MAT[(j, i)] * _totalVars[i] for i in INTVARS) +
                    sum(MAT[(j, i)] * _totalVars[i] for i in CONVARS))
        
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsStrong[i] for i in INTVARS) 
                    + sum(MATFixed[(j, i)] * contVarsStrong[i] for i in CONVARS) 
                    + sum(MATFixed[(j, i)] * slackVarsStrong[i] for i in SLACKVARS)  == RHS[j])
    
    m.optimize()

    status = m.getAttr('Status')
    # print(status)

    if status in [3, 4]:
        return _totalVars
    
    temp = {}
    for k, v in intVarsStrong.items():
        temp[k] = changeValue(v.X)
    for k, v in contVarsStrong.items():
        temp[k] = changeValue(v.X)
    
    return temp


# In[18]:


def RVFSubproblem(currentTheta, timelimit=None):

    m = Model()
    if timelimit:
        m.setParam("TimeLimit", timelimit);
    # m.setParam("LogToConsole", 0);
    m.params.SolutionLimit = 1
    m.params.NonConvex = 2
    
    # set it to True if you want to use presolve option of gurobi
    presolve = False
    
    if presolve:
        m.Params.Presolve = 2

    if presolve:
        dualLB = -GRB.INFINITY
    else:
        dualLB = -1e20

    thetaVar = m.addVar(vtype = GRB.CONTINUOUS, name = "theta")
    intVars = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "int_var")
    conVars = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "cont_var")
    slackVars = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack_var")
    dualVarsVarRHS = m.addVars(len(CONSVARRHS), len(intPartList), vtype = GRB.CONTINUOUS, lb = dualLB, ub = 0,
                               name = "dual_varying_RHS")
    dualVarsFixedRHS = m.addVars(len(CONSFIXEDRHS), len(intPartList), vtype = GRB.CONTINUOUS, lb = dualLB,
                               name = "dual_fixed_RHS")
            
    m.setObjective(1, GRB.MAXIMIZE)

    for k in range(len(intPartList)):
        m.addConstr(thetaVar <=
              sum(OBJ[i]*(intPartList[k][i] - intVars[i]) for i in INTVARS)  
            - sum(OBJ[j]*conVars[j] for j in CONVARS)
            + sum(MAT[(i, j)]*dualVarsVarRHS[(i, k)]*(intVars[j] - intPartList[k][j]) for j in INTVARS for i in CONSVARRHS)
            + sum(MAT[(i, j)]*dualVarsVarRHS[(i, k)]*(conVars[j]) for j in CONVARS for i in CONSVARRHS)
            - sum(MATFixed[(i, j)]*dualVarsFixedRHS[(i, k)]*(intPartList[k][j]) for j in INTVARS for i in CONSFIXEDRHS)
            + sum(RHS[i]*dualVarsFixedRHS[(i, k)] for i in CONSFIXEDRHS))

    m.addConstr(thetaVar + sum(OBJ[i]*intVars[i] for i in INTVARS) + sum(OBJ[i]*conVars[i] for i in CONVARS) <= U + 1)

    for k in range(len(intPartList)):
        for j in CONVARS:
            m.addConstr(sum(MAT[(i, j)]*dualVarsVarRHS[(i, k)] for i in CONSVARRHS)
                          + sum(MATFixed[(i, j)]*dualVarsFixedRHS[(i, k)] for i in CONSFIXEDRHS) <= OBJ[j])
    
    for k in range(len(intPartList)):
        for j in SLACKVARS:
            m.addConstr(sum(MATFixed[(i, j)]*dualVarsFixedRHS[(i, k)] for i in CONSFIXEDRHS) <= 0)

    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVars[i] for i in INTVARS)
                      + sum(MATFixed[(j, i)] * conVars[i] for i in CONVARS)
                      + sum(MATFixed[(j, i)] * slackVars[i] for i in SLACKVARS) == RHS[j])

    
    # print('currentTheta', currentTheta)
    m.addConstr(thetaVar >= currentTheta, name = "theta_threshold")

    m.optimize() 
    status = m.getAttr('Status')
    
    if status in [3, 4]:  
        # print('in CR: infeasible')
        return -1, currentTheta
    
    if status == 9:
        return -2, currentTheta
     
    total_part = convertWeakToStrongNDP({**dict((k, round(changeValue(v.X))) for k,v in intVars.items()), **dict((k, round(changeValue(v.X))) for k,v in conVars.items())}, _print=False)
   
    temp_ndp = () 
    temp_ndp = temp_ndp + (round(sum(OBJ[j]*total_part[j] for j in INTVARS)) +
                           round(sum(OBJ[j]*total_part[j] for j in CONVARS), 5),)
    
    for k in CONSVARRHS:
        temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*total_part[l] for l in INTVARS)) +
                               round(sum(MAT[(k, l)]*total_part[l] for l in CONVARS),5),) 

    # print('intPartList', intPartList)
    # print('total_part', {i: total_part[i] for i in INTVARS})
    if {i: total_part[i] for i in INTVARS} in intPartList:
        # print('int part already exists!')
        currentTheta = currentTheta * 10
        # print('currentTheta inside', currentTheta)
        theta_threshold_c = m.getConstrByName('theta_threshold')
        m.remove(theta_threshold_c)
        m.addConstr(thetaVar >= currentTheta, name = "theta_threshold")

    else:
        # print('new int part')
        currentTheta = defaultTheta
        EF.append(temp_ndp)    
        intPartList.append({k: round(total_part[k]) for k in INTVARS})
        contPartList.append({k:total_part[k] for k in CONVARS})
    # print('in CR, feasibility status {}'.format(status))
    return None, currentTheta


# In[19]:


# Generate CRs for each int parts
def generateCR(_int, _sampleZeta):
    m = Model()
    m.setParam("LogToConsole", 0);

    contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(OBJ[i] * contVarsInit[i] for i in CONVARS), GRB.MINIMIZE)

    RHS_int = []
    for j in CONSFIXEDRHS:
        RHS_int.append(sum(MATFixed[(j, i)] * _int[i] for i in INTVARS))
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j] - RHS_int[j])

    RHS_int_var = []
    for j in CONSVARRHS:
        RHS_int_var.append(sum(MAT[(j, i)] * _int[i] for i in INTVARS))
    
    for j in CONSVARRHS:
        m.addConstr(sum(MAT[(j, i)] * contVarsInit[i] for i in CONVARS) <= _sampleZeta[j] - RHS_int_var[j])
        
    m.optimize()
    
    status = m.getAttr('Status')
    if status in [3,4]:
        return U
        
    objValInt = sum(OBJ[i] * _int[i] for i in INTVARS) 
    
    return (m.objVal + objValInt)


# In[20]:


def calcRVFUBAtZeta(sampleZeta):
    UB_temp_lst = []
    for int_parts in intPartList:
        UB_temp_lst.append(generateCR(int_parts, sampleZeta))
        
    return min(UB_temp_lst)


# In[21]:


def findGapOfZeta(sampleZeta):
    
    UB_bound = calcRVFUBAtZeta(sampleZeta)

    st = time.time()

    LB_bound = forest.evaluate_forest_dual_function(sampleZeta)

    et = time.time()

    elapsed = et - st
    global totTimeDfEval

    totTimeDfEval += elapsed

    if LB_bound >= 1e+20:
        return UB_bound, None
            
    return UB_bound, LB_bound


# In[22]:


def generate_sequence(lst):
    result = []
    result.extend([lst[0], lst[-1]])
    # print('result', result)
    result.append((lst[0] + lst[-1]) / 2)
    result = sorted(result)

    for i in range(2):
        temp_lst = []
        for j in range(len(result)-1):
            mid_value = (result[j] + result[j+1]) / 2
            temp_lst.append(mid_value)

        result.extend(temp_lst)
        
        # while len(temp_lst) >= 2:
        #     result.append(temp_lst.pop(0))
        #     result.append(temp_lst.pop(-1))
            
        # if len(temp_lst) == 1:
        #     result.append(temp_lst.pop(0))
        
        result = sorted(result)
        
    return result


# In[23]:


Zeta_LB = []
for i in CONSVARRHS:
    Zeta_LB.append(generateLBZeta(i))

Zeta_UB = []
for i in CONSVARRHS:
    Zeta_UB.append(generateUBZeta(i))

tempListZeta = []

for i in CONSVARRHS:
    my_list = []
    my_list.append(Zeta_LB[i])
    my_list.append(Zeta_UB[i])
    tempListZeta.append(generate_sequence(my_list))

tempListZeta = list(itertools.product(*tempListZeta))

print(Zeta_LB)
print(Zeta_UB)

tempListZeta = generate_samples_hypercube(num_samples, Zeta_LB, Zeta_UB, len(CONSVARRHS), method='lhs')

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
contPartList = []
EF = []
mainIter = 0
RVFIter = 0 

iterLimit = len(tempListZeta)

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
    total_part = convertWeakToStrongNDP(opt_sol, _print=False)
    temp_ndp = ()
    temp_ndp = temp_ndp + (round(sum(OBJ[j]*total_part[j] for j in INTVARS)) + round(sum(OBJ[j]*total_part[j] for j in CONVARS)),)

    for k in CONSVARRHS: 
        temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*total_part[l] for l in INTVARS)) + round(sum(MAT[(k, l)]*total_part[l] for l in CONVARS)),)

    if temp_ndp not in EF:
        # print('New NDP (SYMPHONY): ', temp_ndp)
        EF.append(temp_ndp)
        intPartList.append({i:round(total_part[i]) for i in INTVARS})
        contPartList.append({i:total_part[i] for i in CONVARS})


feas_status = 0
EF_len_main = 0
EF_len_after_main = 0
numLBGRVF = 0
LBGRVFLst = []
timeSurpassed = False

for sample_zeta in tempListZeta:
    print('------------------------')
    # print('len(EF)', len(EF))
    print('totalIter', mainIter)
    # print('iterLimit', iterLimit)
    # Convert to a python list 
    sample_zeta = sample_zeta.astype(float).tolist()
    mainIter += 1
    
    UB, LB = findGapOfZeta(sample_zeta)

    if LB == None:
        continue

    gap_zeta = UB - LB

    if (gap_zeta <= 0.1):
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
        total_part = convertWeakToStrongNDP(opt_sol, _print=False)
        temp_ndp = ()
        temp_ndp = temp_ndp + (round(sum(OBJ[j]*total_part[j] for j in INTVARS)) + round(sum(OBJ[j]*total_part[j] for j in CONVARS)),)
    
        for k in CONSVARRHS: 
            temp_ndp = temp_ndp + (round(sum(MAT[(k, l)]*total_part[l] for l in INTVARS)) + round(sum(MAT[(k, l)]*total_part[l] for l in CONVARS)),)

        if {i: total_part[i] for i in INTVARS} not in intPartList:
            # print('New NDP (SYMPHONY): ', temp_ndp)
            EF.append(temp_ndp)
            intPartList.append({i:round(total_part[i]) for i in INTVARS})
            contPartList.append({i:total_part[i] for i in CONVARS})

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
    feas_status, currentTheta = RVFSubproblem(currentTheta, timelimit=remainingTime)
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
print('                 Continuous: ', numContVars)
print('Total number of constraints: ', numConsFixed)
print(' Total number of objectives: ', numConsVar + 1)

print('---------------------------------------------')
print('Statistics')
print('---------------------------------------------')

print('Total Running Time (sec): ', totalTime)
print('Total Running Time w/o LPs (sec): ', totalTimeWOLP)
print('EF: ' , EF)
print('Number of total stability regions: ', len(EF))
print('Number of main algo iterations: ', mainIter)
print('Number of iterations that RVF called: ', RVFIter)
print('Number of total iterations: ', mainIter + RVFIter)
print('Len EF inside the main: ', EF_len_main)
print('Len EF after the main: ', EF_len_after_main)
print('Number of cases that LB is greater than RVF: ', numLBGRVF)
print('List of zetas that LB is greater than RVF: ', LBGRVFLst)
print('Number of Total Nodes of the BB: ', sym.get_tree_size())
print('---------------------------------------------')
print('Forest Statistics')
print('---------------------------------------------')
print('Total number of trees: ', len(forest.trees))
print('Total time evaluating dual functions (sec): ', totTimeDfEval)
print('Total time selecting trees (sec): ', totTimeSelectTree)

# with open("Results_MILP_multiobj_details_{}.txt".format(instance.split('.')[0]), "a") as _filed:
#     if timeSurpassed:
#         _filed.write('Finished due to time limit!' + '\n') 
#     else: 
#         _filed.write('Finished!' + '\n') 
        
#     _filed.write('Total Running Time (sec): ' + str(totalTime) + '\n') 
#     _filed.write('Total Running Time w/o LPs (sec): ' + str(totalTimeWOLP) + '\n') 

#     if timeSurpassed:
#         _filed.write('Approximate EF: ' + str(EF) + '\n') 
#     else: 
#         _filed.write('EF: ' + str(EF) + '\n') 
        
#     _filed.write('Number of total stability regions: ' + str(len(EF)) + '\n') 
#     _filed.write('Number of main algo iterations: ' + str(mainIter) + '\n') 
#     _filed.write('Number of iterations that RVF called: ' + str(RVFIter) + '\n') 
#     _filed.write('Number of total iterations: ' + str(mainIter + RVFIter) + '\n') 
#     _filed.write('Len EF inside the main: ' + str(EF_len_main) + '\n') 
#     _filed.write('Len EF after the main: ' + str(EF_len_after_main) + '\n') 
#     _filed.write('Number of cases that LB is greater than RVF: ' + str(numLBGRVF) + '\n') 
#     _filed.write('List of zetas that LB is greater than RVF: ' + str(LBGRVFLst) + '\n') 
#     _filed.write('Number of Total Nodes of the BB: ' + str(sym.get_tree_size())) 


# In[ ]:


del sym




