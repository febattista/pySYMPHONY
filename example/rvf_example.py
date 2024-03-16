import os, sys

sys.path.append("..")

import numpy as np
from symphony import *

# Symphony's Shared library absolute path
lib_path = "/home/feb223/rvf/build-SYM-debug/lib/libSym.so"

Symphony.dlopen(lib_path)

# MPS file
model_path = "./RVF_example.MPS"

# RHS vector (m = 1)
zeta_lst = np.linspace(-56.5, 5, 2000)

# Create SYMPHONY environments
sym_df = Symphony()
sym_rvf = Symphony()

# Set additional parameters
sym_df.set_param("verbosity", -2)
sym_rvf.set_param("verbosity", -2)

# Load the problem
if (sym_df.read_mps(model_path)== FUNCTION_TERMINATED_ABNORMALLY):
    print("Unable to read the MPS file!")
    exit(1)

if (sym_rvf.read_mps(model_path)== FUNCTION_TERMINATED_ABNORMALLY):
    print("Unable to read the MPS file!")
    exit(1)

# Enable Warm Start
sym_df.enable_warm_start()
sym_rvf.enable_warm_start()

# Compute the Full RVF
print("========================")
print("     Computing RVF      ")
print("========================")
print("Building the RVF for zeta = %.2f, ..., %.2f" 
      % (zeta_lst[0], zeta_lst[-1]))
print("Solving...")

rvf_val_lst = []
for zeta in zeta_lst:
    sym_rvf.set_row_upper(0, zeta)
    if (sym_rvf.warm_solve() == FUNCTION_TERMINATED_ABNORMALLY):
        print("Error solving RVF!")
        exit(1)
    # Succesfully solved
    objval = sym_rvf.get_obj_val()
    rvf_val_lst.append(objval)

print("Done!")


# Build the dual function by evaluating the RVF
# at a bunch of RHSs
print("========================")
print("   Building Dual Func   ")
print("========================")
rhss = [40/9, -6, -11, -40]

print("Building a DF for zeta =", rhss) 
print("Solving...")

for rhs in rhss:
    sym_df.set_row_upper(0, rhs)
    if (sym_df.warm_solve() == FUNCTION_TERMINATED_ABNORMALLY):
        print("Error solving RVF!")
        exit(1)
    # Succesfully solved
    # Create and update the dual function
    sym_df.build_dual_function()

print("Done!")

print("========================")
print("  Evaluating Dual Func  ")
print("========================")
print("Evaluating the Dual Func at zeta = %.2f, ..., %.2f" 
      % (zeta_lst[0], zeta_lst[-1]))

df_val_lst = []

for zeta in zeta_lst:
    # Evaluate the dual function
    dual = sym_df.evaluate_dual_function([zeta])
    df_val_lst.append(dual)

print("Done!")

# Assert the dual function is strong at every RHS
for i, _ in enumerate(zeta_lst):
    assert abs(df_val_lst[i] - rvf_val_lst[i]) <= 1e-6

print("Dual function is strong at all RHS!")

# Close the environments 
del sym_df, sym_rvf
