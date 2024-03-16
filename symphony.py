import os
import numpy as np
from cffi import FFI

#-------------------- C functions and struct declarations ---------------------#
cdefs = """ 
    typedef struct MIPDESC MIPdesc;
    typedef struct WARM_START_DESC warm_start_desc;
    typedef struct SYM_ENVIRONMENT sym_environment;
    typedef struct CUT_POOL cut_pool;

    void sym_version(void);
    
    sym_environment * sym_open_environment(void);
    
    int sym_close_environment(sym_environment *env);
    
    int sym_reset_environment(sym_environment *env);
    
    int sym_set_defaults(sym_environment *env);
    
    int sym_parse_command_line(sym_environment *env, int argc, char **argv);
    
    int sym_set_user_data(sym_environment *env, void *user);
    
    int sym_get_user_data(sym_environment *env, void **user);
    
    int sym_read_mps(sym_environment *env, char *infile);
    
    int sym_read_lp(sym_environment *env, char *infile);

    int sym_read_gmpl(sym_environment *env, char *modelfile, 
                                        char *datafile);
    
    int sym_write_mps(sym_environment *env, char *infile);
    
    int sym_write_lp(sym_environment *env, char *infile);
    
    int sym_load_problem(sym_environment *env);
    
    int sym_find_initial_bounds(sym_environment *env);
    
    int sym_solve(sym_environment *env);
    
    int sym_warm_solve(sym_environment *env);
    
    int sym_mc_solve(sym_environment *env);
    
    int sym_create_permanent_cut_pools(sym_environment *env,
                                                int *cp_num); 
    
    cut_pool **sym_get_permanent_cut_pools(sym_environment *env);
    
    int sym_explicit_load_problem(sym_environment
                            *env, int numcols, 
                    int numrows, int *start, int *index, 
                    double *value, double *collb,
                    double *colub, char *is_int, double *obj,
                    double *obj2, char *rowsen,
                    double *rowrhs, double *rowrng,
                    char make_copy);   

    int sym_is_abandoned(sym_environment *env);

    int sym_is_proven_optimal(sym_environment *env);
    
    int sym_is_proven_primal_infeasible(sym_environment *env);	  
    
    int sym_is_iteration_limit_reached(sym_environment *env); 
    
    int sym_is_time_limit_reached(sym_environment *env);
    
    int sym_is_target_gap_achieved(sym_environment *env);

    int sym_get_status(sym_environment *env);
    
    int sym_get_num_cols(sym_environment *env, int *numcols);
    
    int sym_get_num_rows(sym_environment *env, int *numrows);
    
    int sym_get_num_elements(sym_environment *env,
                                        int *numelems);
    
    int sym_get_col_lower(sym_environment *env,
                                        double *collb);
    
    int sym_get_col_upper(sym_environment *env,
                                        double *colub);
    
    int sym_get_row_sense(sym_environment *env,
                                        char *rowsen);
    
    int sym_get_rhs(sym_environment *env,
                                        double *rowrhs);

    int sym_get_matrix(sym_environment *env, int *nz,
                                        int *matbeg, 
                                        int *matind, double *matval);
    
    int sym_get_row_range(sym_environment *env,
                                            double *rowrng);
    
    int sym_get_row_lower(sym_environment *env,
                                            double *rowlb);
    
    int sym_get_row_upper(sym_environment *env,
                                            double *rowub);
    
    int sym_get_obj_coeff(sym_environment *env,
                                            double *obj);
    
    int sym_get_obj2_coeff(sym_environment *env,
                                            double *obj2);
    
    int sym_get_obj_sense(sym_environment *env, int *sense);
    
    int sym_is_continuous(sym_environment *env, int index,
                                            int *value);
    
    int sym_is_binary(sym_environment *env, int index,
                                        int *value);
    
    int sym_is_integer(sym_environment *env, int index,
                                        char *value);
    
    double sym_get_infinity();

    int sym_get_col_solution(sym_environment *env,
                                        double *colsol);
    
    int sym_get_sp_size(sym_environment *env, int *size);
    
    int sym_get_sp_solution(sym_environment *env, int index,
                                        double *colsol, double *objval);
    
    int sym_get_row_activity(sym_environment *env, double *rowact);
    
    int sym_get_obj_val(sym_environment *env,
                                        double *objval);
    
    int sym_get_primal_bound(sym_environment *env,
                                                double *ub);
    
    int sym_get_iteration_count(sym_environment *env,
                                                int *numnodes);
    
    int sym_set_obj_coeff(sym_environment *env,
                                            int index, double value);
    
    int sym_set_obj2_coeff(sym_environment *env,
                                            int index, double value);
    
    int sym_set_col_lower(sym_environment *env,
                                            int index, double value);
    
    int sym_set_col_upper(sym_environment *env,
                                            int index, double value);
    
    int sym_set_row_lower(sym_environment *env,
                                            int index, double value);
    
    int sym_set_row_upper(sym_environment *env,
                                            int index, double value);
    
    int sym_set_row_type(sym_environment *env,
                                        int index, char rowsense, 
                                        double rowrhs, double rowrng);
    
    int sym_set_obj_sense(sym_environment *env, int sense);
    
    int sym_set_col_solution(sym_environment *env,
                                            double * colsol);
    
    int sym_set_primal_bound(sym_environment *env,
                                            double bound);
    
    int sym_set_continuous(sym_environment *env, int index);
    
    int sym_set_integer(sym_environment *env, int index);
    
    int sym_set_col_names(sym_environment *env,
                                            char **colname);
    
    int sym_add_col(sym_environment *env, int numelems,
                                    int *indices, double *elements,
                                    double collb, double colub,
                                    double obj, char is_int, char *name);
    
    int sym_add_row(sym_environment *env, int numelems,
                                    int *indices, double *elements,
                                    char rowsen, double rowrhs,
                                    double rowrng);
    
    int sym_delete_cols(sym_environment *env, int num,
                                        int * indices);
    
    int sym_delete_rows(sym_environment *env, int num,
                                        int * indices);
    
    int sym_write_warm_start_desc(warm_start_desc *ws,
                                            char *file);
    
    warm_start_desc * sym_read_warm_start(char *file);
    
    void sym_delete_warm_start(warm_start_desc *ws);
    
    warm_start_desc * sym_get_warm_start(sym_environment *env, 
                                            int copy_warm_start);
    
    int sym_set_warm_start(sym_environment *env,
                                            warm_start_desc *ws);
    
    int sym_set_int_param(sym_environment *env,
                                            const char *key, int value);
    
    int sym_set_dbl_param(sym_environment *env,
                                            const char *key, double value);
    
    int sym_set_str_param(sym_environment *env,
                                            const char *key,
                                            const char *value);
    
    int sym_get_int_param(sym_environment *env,
                                            const char *key, int *value);
    
    int sym_get_dbl_param(sym_environment *env,
                                            const char *key, double *value);
    
    int sym_get_str_param(sym_environment *env,
                                            const char *key, char **value);
    
    int sym_get_lb_for_new_rhs(sym_environment *env,
                                            int rhs_cnt, int *new_rhs_ind,
                                            double *new_rhs_val,
                                            int lb_cnt, int *new_lb_ind,
                                            double *new_lb_val,
                                            int ub_cnt, int *new_ub_ind,
                                            double *new_ub_val,
                                            double *lb_for_new_rhs);
    
    int sym_get_dual_pruned(sym_environment *env,
                                        double ** dual_pieces,
                                        int* num_pieces,
                                        int MAX_ALLOWABLE_NUM_PIECES);
    
    int sym_get_ub_for_new_rhs(sym_environment *env,
                                        int cnt, int *new_rhs_ind,
                                        double *new_rhs_val,
                                        double *ub_for_new_rhs);

    int sym_get_ub_for_new_obj(sym_environment *env,
                                                int cnt, int *new_obj_ind,
                                                double *new_obj_val,
                                                double *ub_for_new_obj);
    
    warm_start_desc * sym_create_copy_warm_start(warm_start_desc * ws); 
    
    MIPdesc * sym_create_copy_mip_desc(sym_environment *env);
    
    MIPdesc * sym_get_presolved_mip_desc(sym_environment *env); 
    
    sym_environment * sym_create_copy_environment(sym_environment *env);
    
    int sym_test(sym_environment *env, int argc,
                                    char **argv, int *test_status);
    
    void sym_print_statistics(sym_environment *env,
                                            double start_time,
                                            double finish_time);
    
    double sym_wall_clock(double *T);
    
    int sym_set_param(sym_environment *env, char *line);
    
    int sym_free_env(sym_environment *env);

    void print_dual_function(sym_environment *env);

    int sym_build_dual_func(sym_environment *env);

    int sym_evaluate_dual_function(sym_environment *env, 
            double *new_rhs, int size_new_rhs, double *dual_bound);
 
    """

#-------------------------------- CONSTANTS -----------------------------------#

FUNCTION_TERMINATED_NORMALLY    = 0
FUNCTION_TERMINATED_ABNORMALLY  = -1

INF = float("inf")

#------------------------------------------------------------------------------#

class Symphony():

    is_dl_loaded = False
    ffi = None
    symlib = None

    # Absolute path to the libSym, extention change based on the OS:
    # - Windows: .dll
    # - Linux: .so
    # - OSX: .dylib
    @staticmethod
    def dlopen(lib_path):
        if not Symphony.is_dl_loaded:
            Symphony.ffi = FFI()
            # Load the shared library
            try:
                Symphony.symlib = Symphony.ffi.dlopen(lib_path)
                Symphony.ffi.cdef(cdefs)
                Symphony.is_dl_loaded = True
            except Exception as e:
                print(f"An exception occurred: {e}")
                print(f"Shared library cannot be loaded.")
                exit(1)
        else:
            print("Shared Library already loaded.")

    def __init__(self):
        if not Symphony.is_dl_loaded:
            raise Exception("Shared Library is not loaded yet.")
        
        self._env = Symphony.symlib.sym_open_environment()
        self.model_file = ""                        # Model file
        self.n = -1                                 # Num cols
        self.m = -1                                 # Num rows
        self.warm_start_is_on = False

    def __del__(self):
        if Symphony.is_dl_loaded:
            if self._env:
                return Symphony.symlib.sym_close_environment(self._env)

    def read_mps(self, model: str):
        termcode = FUNCTION_TERMINATED_ABNORMALLY
        # Check if the path exists and ends with ".mps"
        if os.path.exists(model) and model.lower().endswith(".mps"):
            arg_file = Symphony.ffi.new("char []", str.encode(model))
            termcode = Symphony.symlib.sym_read_mps(self._env, arg_file)
            self.model_file = model
            var = Symphony.ffi.new("int *")
            Symphony.symlib.sym_get_num_cols(self._env, var)
            self.n = var[0]
            Symphony.symlib.sym_get_num_rows(self._env, var)
            self.m = var[0]
        return termcode
    
    def write_lp(self, filename):
        arg_file = Symphony.ffi.new("char []", str.encode(filename))
        termcode = Symphony.symlib.sym_write_lp(self._env, arg_file)
        return termcode

    def set_param(self, key: str, value):
        termcode = FUNCTION_TERMINATED_ABNORMALLY
        # Int param
        if type(value) == int:
            termcode = Symphony.symlib.sym_set_int_param(self._env, 
                                                str.encode(key), value)
        # Bool param
        elif type(value) == bool:
            termcode = Symphony.symlib.sym_set_int_param(self._env, 
                                            str.encode(key), 1 if value else 0)
        # Float param
        elif type(value) == float:
            termcode = Symphony.symlib.sym_set_dbl_param(self._env, 
                                                str.encode(key), value)
        # String param
        elif type(value) == str:
            termcode = Symphony.symlib.sym_set_str_param(self._env, str.encode(key), 
                                        str.encode(value))
        else:
            print("Not a valid parameter value.")
        return termcode

    def enable_warm_start(self):
        # Set parameters needed for warm start to work
        if not self.warm_start_is_on:
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("keep_warm_start"), True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("keep_dual_function_description"), 
                                     True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("should_use_rel_br"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("use_hot_starts"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("should_warmstart_node"), True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("sensitivity_analysis"), True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("sensitivity_bounds"), True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("sensitivity_rhs"), True)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("set_obj_upper_lim"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("do_primal_heuristic"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("prep_level"), -1)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("tighten_root_bounds"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("max_sp_size"), 100)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("do_reduced_cost_fixing"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("generate_cgl_cuts"), False)
            Symphony.symlib.sym_set_int_param(self._env, 
                                     str.encode("max_active_nodes"), 1)
            self.warm_start_is_on = True
        else:
            print("Warm start is already enabled.")

    def solve(self):
        termcode = Symphony.symlib.sym_solve(self._env)
        return termcode
    
    def warm_solve(self):
        if self.warm_start_is_on:
            termcode = Symphony.symlib.sym_warm_solve(self._env)
            return termcode
        else:
            print("Warm start is not enabled.")
            return FUNCTION_TERMINATED_ABNORMALLY

    def is_proven_optimal(self):
        return Symphony.symlib.sym_is_proven_optimal(self._env)

    def is_proven_primal_infeasible(self):
        return Symphony.symlib.sym_is_proven_primal_infeasible(self._env)

        
    def get_obj_val(self):
        objval = Symphony.ffi.new("double *")
        termcode = Symphony.symlib.sym_get_obj_val(self._env, objval)
        if termcode == FUNCTION_TERMINATED_ABNORMALLY:
            return -INF
        else:
            return objval[0]
    
    def get_col_solution(self):
        colsol = Symphony.ffi.new("double[%d]" % self.n)
        termcode = Symphony.symlib.sym_get_col_solution(self._env, colsol)
        if termcode == FUNCTION_TERMINATED_ABNORMALLY:
            return []
        else:
            return list(colsol)

    def set_row_lower(self, index: int, rhs: float):
        return Symphony.symlib.sym_set_row_lower(self._env, index, rhs)

    def set_row_upper(self, index: int, rhs: float):
        return Symphony.symlib.sym_set_row_upper(self._env, index, rhs)
    
    def print_dual_function(self):
        Symphony.symlib.print_dual_function(self._env)
    
    def build_dual_function(self):
        termcode = Symphony.symlib.sym_build_dual_func(self._env)
        return termcode
    
    def evaluate_dual_function(self, new_rhs):
        dual_bound = Symphony.ffi.new("double *")
        termcode = Symphony.symlib.sym_evaluate_dual_function(self._env, new_rhs, len(new_rhs), dual_bound)
        if termcode == FUNCTION_TERMINATED_ABNORMALLY:
            return -INF
        else:
            return dual_bound[0]

