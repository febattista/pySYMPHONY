import sys, itertools, time
from symphony import *
from gurobipy import *


# This is the path pointing to the dynamic library
lib_path = r'/home/feb223/rvf/build-SYMPHONY-rvf/lib/libSym.so'

INTVARS      = []
CONVARS      = []
SLACKVARS    = []
MAT          = []
MATFixed     = []
RHS          = []
CONSFIXEDRHS = []
CONSVARRHS   = []
UB_I         = []


def parse_command_line_args():
    args_dict = {}
    i = 1  # Start from the second argument since the first one is the script name
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg.startswith("-"):
            if i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith("-"):
                key = arg[1:]  # Remove the leading '-' character
                value = sys.argv[i + 1]
                args_dict[key] = value
                i += 1  # Skip the next argument as it's the value
            else:
                print(f"Missing value for option {arg}.")
        else:
            print(f"Illegal argument format: {arg}")
        i += 1
    return args_dict

# Read KP data file
def readKPdat(filedat):
    global INTVARS      
    global CONVARS      
    global SLACKVARS   
    global MAT          
    global MATFixed     
    global RHS          
    global CONSFIXEDRHS 
    global CONSVARRHS   
    global UB_I   
    f = open(filedat[:-4] + '.dat', "r")
    content = f.read()   
    cont = content.split("\n")
    cont = [i.rstrip() for i in cont]
    numObj = int(cont[0])
    numVar = int(cont[1])
    capacity = int(cont[2])

    mapping = [(i, ' ') for i in ['[',']',',']]

    MAT = {}
    OBJ = []
    for i in range(numObj):
        tmp = cont[3 + i]
        for k, v in mapping:
            tmp = tmp.replace(k, v)
        tmp = list(map(int, tmp.split())) 
        tmp.append(0)
        if i == 0:
            OBJ = [-tmp[k] for k in range(numVar+1)]
        else:
            MAT.update({(i-1, j):-tmp[j] for j in range(numVar+1)})

    numConst = 1

    MATFixed = {}
    for i in range(numConst):
        tmp = cont[i + 3 + numObj]
        for k, v in mapping:
            tmp = tmp.replace(k, v)
        tmp = list(map(int, tmp.split())) 
        for j in range(numVar):
            MATFixed[(i, j)] = tmp[j] 
    MATFixed[(0, numVar)] = 1  
    RHSList = [capacity]

    # num_NDP = int(cont[1 + 3 + numObj].split('num_NDP =')[1])

    numVars = numVar + numConst
    numIntVars = numVar    
    numConsVar = numObj - 1
    numConsFixed = 1
    INTVARS = range(numIntVars)
    SLACKVARS = range(numIntVars, numVars)
    CONSVARRHS = range(numConsVar)
    CONSFIXEDRHS = range(numConsFixed)

    RHS = {i : RHSList[i] for i in range(0, len(RHSList))}
    UB_I = 1

    fixed_RHS_list = list(RHS.values())

    eps = 0.5

    M = {}

    for i in CONSVARRHS:
        M[i] = sum(-MAT[(i, j)] for j in range(numVar)) + 1

    RHS_list = []
    for i in range(numObj-1):
        RHS_list.append(0) 

    for j in range(len(fixed_RHS_list)):
        RHS_list.append(fixed_RHS_list[j])      

# Read SPP data file
def readSPPdat(filedat):
    global INTVARS      
    global CONVARS      
    global SLACKVARS   
    global MAT          
    global MATFixed     
    global RHS          
    global CONSFIXEDRHS 
    global CONSVARRHS   
    global UB_I         
    f = open(filedat[:-4] + '.dat', "r")
    content = f.read()   
    cont = content.split("\n")
    cont = [i.rstrip() for i in cont]
    # SPP
    numObj = 2
    numConst, numVar = int(cont[0].split()[0]), int(cont[0].split()[1])
    c1_coefs = [int(i) for i in cont[1].split()]
    c2_coefs = [int(i) for i in cont[2].split()]
    cont = cont[3:-1]
    listIndexes = [[int(j) - 1 for j in i.split()] for i in cont[1::2] if i != ''] 
    c1_coefs = [-i for i in c1_coefs] + [0 for i in range(numConst)]
    c2_coefs = [-i for i in c2_coefs] + [0 for i in range(numConst)]

    numVars = numVar + numConst
    numIntVars = numVar 
    numContVars = 0
    numConsVar = 1
    numConsFixed = numConst
    INTVARS = range(numIntVars)
    SLACKVARS = range(numIntVars, numVars)
    VARS = range(numVars)
    CONSVARRHS = range(numConsVar)
    CONSFIXEDRHS = range(numConsFixed)

    OBJ = c1_coefs
    MAT = {(0, i):c2_coefs[i] for i in range(numIntVars)}
    RHS = {i:1 for i in range(numConst)}

    fixed_RHS_list = list(RHS.values())

    MATFixed = {}
    for j in range(numConst):
        for i in range(numVar + numConst):
            MATFixed[(j, i)] = 0

    k = 0
    for j in range(numConst):
        for i in listIndexes[j]:
            MATFixed[(j, i)] = 1
        listIndexes[j].append(numVar + k)
        MATFixed[(j, numVar + k)] = 1
        k += 1
    UB_I = 1

    eps = 0.5

    M = {}

    for i in CONSVARRHS:
        M[i] = sum(-MAT[(i, j)] for j in range(numVar)) + 1

    RHS_list = []
    RHS_list.append(0) 

    for j in range(len(fixed_RHS_list)):
        RHS_list.append(fixed_RHS_list[j])


def generateLBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    if CONVARS:
        contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS) +
                   sum(MAT[(objInd, i)] * contVarsInit[i] for i in CONVARS), GRB.MINIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()

    return m.objVal


def generateUBZeta(objInd):
    m = Model()
    m.setParam("LogToConsole", 0);
    
    intVarsInit = m.addVars(list(INTVARS), vtype = GRB.INTEGER, lb = 0, ub = UB_I, name = "integer variable")
    if CONVARS:
        contVarsInit = m.addVars(list(CONVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "continuous variable")
    slackVarsInit = m.addVars(list(SLACKVARS), vtype = GRB.CONTINUOUS, lb = 0, name = "slack variable")
    
    m.setObjective(sum(MAT[(objInd, i)] * intVarsInit[i] for i in INTVARS) +
                   sum(MAT[(objInd, i)] * contVarsInit[i] for i in CONVARS), GRB.MAXIMIZE)
    
    for j in CONSFIXEDRHS:
        m.addConstr(sum(MATFixed[(j, i)] * intVarsInit[i] for i in INTVARS) +
                    sum(MATFixed[(j, i)] * contVarsInit[i] for i in CONVARS) +
                    sum(MATFixed[(j, i)] * slackVarsInit[i] for i in SLACKVARS) == RHS[j])
    m.optimize()

    return m.objVal


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


def generate_zeta_lst():
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
    return tempListZeta


def solve(model_path, tempListZeta):
    # Load the library by calling this static method
    Symphony.dlopen(lib_path)
    # Create SYMPHONY environment
    sym = Symphony()
    # Set additional parameters
    sym.set_param("verbosity", -2)
    # Load the problem
    sym.read_mps(model_path)
    # Enable Warm Start
    sym.enable_warm_start() 
    st = time.time()
    totalTime = 0
    totalTimeBuild = 0
    totalTimeEval = 0
    num_solve = 0
    if (sym.solve() == FUNCTION_TERMINATED_ABNORMALLY):
        print("Something went wrong!")
    # exit(1)   
    num_solve += 1
    st_build = time.time()
    sym.build_dual_function()
    totalTimeBuild += time.time() - st_build
    sym.print_dual_function()
    
    if sym.is_proven_optimal():
        print("OPTIMAL")
        obj = sym.get_obj_val()
        print("RVF: ", obj)
    else:
        print("INFEASIBLE")
    st_eval = time.time()
    df_val = sym.evaluate_dual_function([])
    totalTimeEval += time.time() - st_eval
    print("DF: ", df_val)
    if (abs(obj - df_val) > 1e-5):
        print("Warning: dual function not strong!")
    
    for sample_zeta in tempListZeta:
        print("-------", sample_zeta ,"--------")
        for i in CONSVARRHS:
            sym.set_row_upper(i, sample_zeta[i])
        sym.warm_solve()
        num_solve += 1
        st_build = time.time()
        sym.build_dual_function()
        totalTimeBuild += time.time() - st_build

        sym.print_dual_function()
        print("Current tree size: ", sym.get_tree_size())
        if sym.is_proven_optimal():
            print("OPTIMAL")
            obj = sym.get_obj_val()
            print("RVF: ", obj)
        else:
            print("INFEASIBLE")
        st_df = time.time()
        df_val = sym.evaluate_dual_function(list(sample_zeta))
        totalTimeEval += time.time() - st_df
        print("DF: ", df_val)
        if (sym.is_proven_optimal() and abs(obj - df_val) > 1e-5):
            print("Warning: dual function not strong!")
        elif (sym.is_proven_primal_infeasible() and df_val < 1e+10):
            print("Warning: dual function not strong!")

    # Take time at the end of the loop
    et = time.time()

    # This is the total time elapsed (including time spent in LPs)
    totalTime += et - st 

    # This is the function giving you the time (in seconds) spent up to this moment in LPs
    print('            Time spent on additional LPs (sec): ', sym.get_lp_cpu_time())
    print('                      Total Running Time (sec): ', totalTime)
    print('              Total Running Time w/o LPs (sec): ', totalTime - sym.get_lp_cpu_time())
    print('   Total Running Time on build_dual_func (sec): ', totalTimeBuild)
    print('     Avg Running Time on build_dual_func (sec): ', totalTimeBuild/num_solve)
    print('Total Running Time on evaluate_dual_func (sec): ', totalTimeEval - sym.get_lp_cpu_time())
    print('  Avg Running Time on evaluate_dual_func (sec): ', (totalTimeEval - sym.get_lp_cpu_time())/num_solve)




def main():
    args_dict = parse_command_line_args()
    tempListZeta = []
    print("Parsed arguments:")
    for key, value in args_dict.items():
        print(f"{key}: {value}")

    if 'F' in args_dict:
        if 'spp' in args_dict['F'].split(os.sep)[-1].lower():
            readSPPdat(args_dict['F'])
            tempListZeta = generate_zeta_lst()
        elif 'mis' in args_dict['F'].split(os.sep)[-1].lower():
            readSPPdat(args_dict['F'])
            tempListZeta = generate_zeta_lst()
        elif 'kp' in args_dict['F'].split(os.sep)[-1].lower():
            readKPdat(args_dict['F'])
            tempListZeta = generate_zeta_lst()
            # print(tempListZeta[:13])
        else:
            exit(1)

        solve(args_dict['F'], tempListZeta)
    
    else:
        exit(2)

if __name__ == "__main__":
    main()