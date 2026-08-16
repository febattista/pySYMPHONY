import sys
import os

import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline
from matplotlib import rcParams
rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['text.usetex'] = True

from symphony import *

# TODO
# - get rid of hard-coded paths
# - automatically create df*.log during run or else incorporate directly into SYMPHONY

lib_path = r'/home/ted/Projects/build-SYMPHONY-rvf-debug/lib/libSym.so'
# Load the library by calling this static method
Symphony.dlopen(lib_path)

def plot_disc(x, y, label, offset_end=0, markersize=2 ,color='b', style='-', line_width=2, threshold=1, zorder=1): 
    discontinuities = np.where(np.abs(np.diff(y)) > threshold)[0]
    start = 0
    for discontinuity in discontinuities:
        end = discontinuity + 1
        plt.plot(x[start:end - offset_end], y[start:end - offset_end], style, color=color, linewidth=line_width, label=label if start == 0 else "", zorder=zorder)

        plt.plot(x[end-1], y[end-1], 'o', color=color, markersize=markersize, fillstyle='none', zorder=zorder) 
        plt.plot(x[end], y[end], 'o', color=color, markersize=markersize, fillstyle='full', zorder=zorder)     

        start = end
        
    plt.plot(x[start:], y[start:], style, color=color, linewidth=line_width, zorder=1, label=label if start == 0 else "")


def compute_rvf(zeta_lst):
    rvf_lst = []
    # Create SYMPHONY environment
    sym = Symphony()
    # Set additional parameters
    sym.set_param("verbosity", -2)
    # Load the problem
    sym.read_mps(MPS_file)

    for zeta in zeta_lst:
        retval = ""
        sym.set_row_upper(0, zeta)
        sym.solve()
        if sym.is_proven_optimal():
            obj = sym.get_obj_val()
            rvf_lst.append(obj)
            retval = "%.5f" % obj
        else:
            rvf_lst.append(np.nan)
            retval = "INF"
        # print(("RHS: %.5f VAL: " + retval) % zeta)

    rvf_lst = np.array(rvf_lst)
    return rvf_lst


def compute_ef(rvf_lst):
    # Compute the efficient frontier from the RVF
    # Find points in which the RVF is flat 
    # and replace with NaN
    ef_lst = np.full_like(rvf_lst, np.nan, dtype=np.float64)
    mask = np.insert(np.abs(np.diff(rvf_lst)) > 1e-4, 0, True)
    ef_lst[mask] = rvf_lst[mask]
    return ef_lst


def compute_df(rhss, zeta_lst, logfile=''):
    df_lst = []
    # Create SYMPHONY environment
    sym = Symphony()
    # Set additional parameters
    sym.set_param("verbosity", -2)
    sym.set_param("do_primal_heuristic", 0)
    # Load the problem
    sym.read_mps(MPS_file)
    sym.enable_warm_start() 

    for rhs in rhss:
        sym.set_row_upper(0, rhs)
        sym.warm_solve()
        sym.build_dual_function()
    
    sym.print_dual_function()

    for zeta in zeta_lst:
        val = sym.evaluate_dual_function([zeta])
        df_lst.append(val)

    df_lst = np.array(df_lst)
    return df_lst



'''
Compute nodes stability region from SYMPHONY's output
produced when the RVF branch is built with -DCHECK_DUAL_FUNC.
There should be an easier way to do this.
'''
def compute_regions(logfile):
    best_node = []
    with open(logfile, "r") as f:
        for line in f:
            if line.startswith("Best disj:"):
                parts = line.strip().split()
                disj_val = int(parts[2].rstrip(','))       # value after Best disj:
                best_node.append(disj_val)
    best_node = np.array(best_node)
    regions = []
    if best_node.size > 0:
        last_val = best_node[0]
        regions.append((last_val, 0))   # always record the first
        for i in range(1, len(best_node)):
            if best_node[i] != last_val:
                regions.append((best_node[i], i))
                last_val = best_node[i]
    print(regions)
    return regions


MPS_file = "/home/ted/Projects/pySYMPHONY/example/RVF_example.MPS"
threshold = 2
line_width = 3
filename = "rvf_example.pdf"


zeta_lb = -53
zeta_ub = 3
num_points = 5000

zeta_lst_large = np.linspace(zeta_lb, zeta_ub, num=num_points)

rvf_lst = compute_rvf(zeta_lst_large)
ef_lst = compute_ef(rvf_lst)

df_lst_1 = compute_df([-40/9], zeta_lst_large)
reg_1 = compute_regions('/home/ted/Projects/pySYMPHONY/dflogs-new2/df1.log')
df_lst_2 = compute_df([-40/9, -55.5], zeta_lst_large)
reg_2 = compute_regions('/home/ted/Projects/pySYMPHONY/dflogs-new2/df2.log')
#df_lst_3 = compute_df([-40/9, -55.5, -73/6, -30], zeta_lst_large) # Federico's original example
df_lst_3 = compute_df([-40/9, -55.5, -11], zeta_lst_large)
reg_3 = compute_regions('/home/ted/Projects/pySYMPHONY/dflogs-new2/df3.log')
#df_lst_4 = compute_df([-40/9, -55.5, -73/6, -11, -8, -30], zeta_lst_large) # Federico's original example
df_lst_4 = compute_df([-40/9, -55.5, -11, -5, -30], zeta_lst_large)
reg_4 = compute_regions('/home/ted/Projects/pySYMPHONY/dflogs-new2/df4.log')

#This is for the workaround for creation of df*.log files 
#sys.exit(0)

# ================================================================
#       LARGE PLOTS
# ================================================================

# Base graph RVF + EF
filename = "rvf_example_base_large.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_large, rvf_lst, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

ofs_rvf = 37
mrkr_rvf = 4.5
lw_rvf = 3

ofs_ef = 20
mrkr_ef = 3.5
lw_ef = 1.5

plot_disc(zeta_lst_large, rvf_lst, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_large, ef_lst, r"${\rm EF}$", color='dodgerblue', style='--', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

plt.ylim(-1, 85)
plt.xlim(-55, 1)
# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend()

plt.savefig(filename)
plt.clf()

# RVF + DF1
filename = "rvf_example_large_df_1.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_large, rvf_lst, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_large, rvf_lst, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_large, df_lst_1, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

x_val = -40/9
plt.annotate(
    r"$\zeta_1 = -4$",            # text
    xy=(x_val, 2),             # point to point at
    xytext=(-5, 20),               # text location
    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
    fontsize=12,
    ha="right",
    va="center"
)

for i, r in enumerate(reg_1):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_1) - 1:
        pos = (x_val + zeta_lst_large[reg_1[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)
    plt.text(pos, -2.5, "Node %d" % (r[0]), fontsize=12, ha='center', va='center')

plt.ylim(-5, 85)
plt.xlim(-55, 1)
# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF2
filename = "rvf_example_large_df_2.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_large, rvf_lst, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_large, rvf_lst, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_large, df_lst_2, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

x_val = -50
plt.annotate(
    r"$\zeta_2 = -50$",            # text
    xy=(x_val, 80),             # point to point at
    xytext=(x_val + 5, 20),               # text location
    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
    fontsize=12,
    ha="right",
    va="center"
)

for i, r in enumerate(reg_2):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_2) - 1:
        pos = (x_val + zeta_lst_large[reg_2[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)
    plt.text(pos, -2.5, "Node %d" % (r[0]), fontsize=12, ha='center', va='center')

plt.ylim(-5, 85)
plt.xlim(-55, 1)
# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF3
filename = "rvf_example_large_df_3.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_large, rvf_lst, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_large, rvf_lst, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_large, df_lst_3, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

# x_val = -73/6
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val - 1, 30, r"$\zeta_3 = -\frac{73}{6}$", fontsize=12, ha='right', va='center')

# x_val = -30
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val - 1, 30, r"$\zeta_4 = -30$", fontsize=12, ha='right', va='center')

#x_val = -73/6
#plt.annotate(
#    r"$\zeta_3 = -\frac{73}{6}$",       # text
#    xy=(x_val, 6),             # point to point at
#    xytext=(-20, 2),       # text location
#    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
#    fontsize=12,
#    ha="right",
#    va="center"
#)

x_val = -11
plt.annotate(
    r"$\zeta_5 = -11$",         # text
    xy=(x_val, 7),              # point to point at
    xytext=(-13, 40),           # text location
    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
    fontsize=12,
    ha="right",
    va="center"
)

for i, r in enumerate(reg_3):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_3) - 1:
        pos = (x_val + zeta_lst_large[reg_3[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)
    # The "+ 3" here is because the tree is complete and has 7 nodes. We only want to consider the leaves
    plt.text(pos, -2.5, "Node %d" % (r[0] + 3), fontsize=5, ha='center', va='center')

plt.ylim(-5, 85)
plt.xlim(-55, 1)
# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF4
filename = "rvf_example_large_df_4.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_large, rvf_lst, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_large, rvf_lst, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_large, df_lst_4, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

# x_val = -11
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val - 1, 30, r"$\zeta_5 = -11$", fontsize=12, ha='right', va='center')

# x_val = -8
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val + 1, 30, r"$\zeta_6 = -8$", fontsize=12, ha='left', va='center')

#x_val = -8
#plt.annotate(
#    r"$\zeta_6 = -8$",          # text
#    xy=(x_val, 7),              # point to point at
#    xytext=(0.5, 35),           # text location
#    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
#    fontsize=12,
#    ha="right",
#    va="center"
#)

x_val = -5
plt.annotate(
    r"$\zeta_1 = -5$",            # text
    xy=(x_val, 2),             # point to point at
    xytext=(-5, 20),               # text location
    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
    fontsize=12,
    ha="right",
    va="center"
)

x_val = -30
plt.annotate(
    r"$\zeta_4 = -30$",       # text
    xy=(x_val, 39),             # point to point at
    xytext=(-40, 10),       # text location
    arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
    fontsize=12,
    ha="right",
    va="center"
)

for i, r in enumerate(reg_4):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_4) - 1:
        pos = (x_val + zeta_lst_large[reg_4[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)


plt.text(-30, -2.5, "Node 4", fontsize=12, ha='center', va='center')
plt.text(  0, -2.5, "Node 3", fontsize=12, ha='right', va='center')

plt.ylim(-5, 85)
plt.xlim(-55, 1)
# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()


# ================================================================
#       SMALL PLOTS
# ================================================================
zeta_lb = -20

ofs_rvf = 30
mrkr_rvf = 10
lw_rvf = 4

ofs_ef = 20
mrkr_ef = 8
lw_ef = 2.5


zeta_lst_small = zeta_lst_large[zeta_lst_large >= zeta_lb]
rvf_lst_small = rvf_lst[zeta_lst_large >= zeta_lb]
ef_lst_small = ef_lst[zeta_lst_large >= zeta_lb]

df_lst_1_small = df_lst_1[zeta_lst_large >= zeta_lb]
df_lst_2_small = df_lst_2[zeta_lst_large >= zeta_lb]
df_lst_3_small = df_lst_3[zeta_lst_large >= zeta_lb]
df_lst_4_small = df_lst_4[zeta_lst_large >= zeta_lb]

# Base graph RVF + EF
filename = "rvf_example_base_small.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_small, rvf_lst_small, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_small, rvf_lst_small, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_small, ef_lst_small, r"${\rm EF}$", color='dodgerblue', style='--', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

line_width = 1
plt.ylim(-10, 30)
plt.xlim(-16, 1)

# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF1
filename = "rvf_example_small_df_1.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_small, rvf_lst_small, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_small, rvf_lst_small, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_small, df_lst_1_small, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

x_val = -40/9
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(-5, 20, r"$\zeta_1 = -5$", fontsize=12, ha='right', va='center')

# plt.annotate(
#     r"$\zeta_1 = -5$",            # text
#     xy=(x_val, 2),             # point to point at
#     xytext=(-5, 20),               # text location
#     arrowprops=dict(arrowstyle="->", color="black", alpha=0.7),
#     fontsize=12,
#     ha="right",
#     va="center"
# )

plt.text(-7, -2.5, "Node 0", fontsize=12, ha='center', va='center')

plt.ylim(-4, 30)
plt.xlim(-16, 1)

# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF2
filename = "rvf_example_small_df_2.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_small, rvf_lst_small, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_small, rvf_lst_small, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_small, df_lst_2_small, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

plt.text(-7, -2.5, "Node 0", fontsize=12, ha='center', va='center')

plt.ylim(-4, 30)
plt.xlim(-16, 1)

# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF3
filename = "rvf_example_small_df_3.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_small, rvf_lst_small, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_small, rvf_lst_small, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_small, df_lst_3_small, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

for i, r in enumerate(reg_3):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_3) - 1:
        pos = (x_val + zeta_lst_large[reg_3[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)
    # The "+ 3" here is because the tree is complete and has 7 nodes. We only want to consider the leaves
    plt.text(pos, -2.5, "Node %d" % (r[0] + 3), fontsize=12, ha='center', va='center')

plt.ylim(-4, 30)
plt.xlim(-16, 1)

# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend(loc="upper right")

plt.savefig(filename)
plt.clf()

# RVF + DF4
filename = "rvf_example_small_df_4.pdf"
# Fill the epigraph
plt.fill_between(zeta_lst_small, rvf_lst_small, y2=100, color='#FFDAB9', alpha=1, label="epigraph")

plot_disc(zeta_lst_small, rvf_lst_small, r"${\rm RVF}$", color='black', style='-', line_width=lw_rvf, threshold=threshold, offset_end=ofs_rvf, markersize=mrkr_rvf)
plot_disc(zeta_lst_small, df_lst_4_small, r"$\underline{\phi}$", color='red', style='-', line_width=lw_ef, threshold=threshold, offset_end=ofs_ef, markersize=mrkr_ef)

# x_val = -11
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val - .5, 20, r"$\zeta_5 = -11$", fontsize=12, ha='right', va='center')

# x_val = -8
# plt.axvline(x=x_val, linestyle="dashed", color="black", alpha=0.3)
# plt.text(x_val + .5, 20, r"$\zeta_6 = -8$", fontsize=12, ha='left', va='center')

for i, r in enumerate(reg_4):
    x_val = zeta_lst_large[r[1]]
    if i < len(reg_4) - 1:
        pos = (x_val + zeta_lst_large[reg_4[i + 1][1]]) / 2
    else: 
        pos = (x_val + (-1)) / 2
    plt.axvline(x=x_val, color="blue", alpha=0.3)

plt.text(-14, -2.5, "Node 4", fontsize=12, ha='center', va='center')
plt.text(-11, -2.5, "Node 6", fontsize=12, ha='center', va='center')
plt.text(-9, -2.5, "Node 5", fontsize=12, ha='center', va='center')
plt.text(-7.1, -2.5, "Node 4", fontsize=9, ha='center', va='center')
plt.text(-3, -2.5, "Node 3", fontsize=12, ha='center', va='center')

plt.ylim(-4, 30)
plt.xlim(-16, 1)

# Optional: Add labels and title
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\phi(\zeta)$")

# Add legend
plt.legend()

plt.savefig(filename)
plt.clf()

# m = 4
# n = 7
# lb = [0, 0, 0, 0, 0, 0, 0]
# ub = [1, 1, np.inf, np.inf, np.inf, np.inf, np.inf]
# rhs = [0, 4, 5, 5]


# T = [
#     {
#         'lb' : lb,
#         'ub' : ub,
#         'leaf' : [0]
#     }
# ]

# duals = [
#     {0: -0.1463414634, 1: 0.1626016260, 4: 2.0162601626, 5: 2.8861788618, 7: 5.3414634146, 8: 9.8211382114, 10: 12.5040650407},
#     {0: -2.3529411765, 1: -2.4117647059, 2: -4.5882352941, 4: -2.7647058824, 5: 14.0588235294, 6: 45.2352941176, 8: 17.1764705882},
#     {0: -1.5339805825, 1: -1.4563106796, 2: -0.9029126214, 4: -0.9902912621, 6: 28.4466019417, 8: 14.4466019417, 10: 4.6407766990},
#     {0: -1.0769230769, 1: -0.9230769231, 2: -0.4153846154, 6: 19.0769230769, 7: 1.5692307692, 8: 12.9230769231, 10: 7.2307692308},
#     {0: -1.7521367521, 1: -2.0512820513, 2: -0.8632478632, 4: -1.8034188034, 6: 35.9829059829, 8: 15.8547008547, 9: 2.0427350427},
#     {0: -1.6129032258, 1: -1.9677419355, 4: -1.5806451613, 5: -3.2580645161, 6: 33.8387096774, 8: 15.5483870968, 9: 2.5161290323},
#     {0: -1351817728.2549221516, 1: -1577120682.6307423115, 2: -6083179772.7163801193, 3: -7660300454.0086593628, 4: -2928938408.8856644630, 5: 61507706586.5989456177, 6: 27712263426.2259025574, 7: 1.5692307692, 8: 4506059103.5164070129, 10: 7.2307681205},
#     {0: -0.6875000000, 1: -0.4687500000, 4: 0.8437500000, 6: 11.0937500000, 7: 2.9062500000, 8: 11.6250000000, 10: 9.4375000000},
#     {0: -1.3333333333, 1: -1.2222222222, 4: -0.5555555556, 5: -3.4444444444, 6: 24.3333333333, 8: 13.7777777778, 10: 5.7777777778},
#     {1: -1.0000000000, 4: 1.0000000000, 5: 9.0000000000, 6: 9.0000000000, 7: 10.0000000000, 8: 12.0000000000, 9: 8.0000000000},
#     {4: 2.0000000000, 5: 5.0000000000, 7: 7.0000000000, 8: 10.0000000000, 9: 2.0000000000, 10: 10.0000000000}
# ]

# for t_idx, t in enumerate(T):
#     for d_idx in t['leaf']:
#         intcpt = 0
#         zeta = duals[d_idx][0] if 0 in duals[d_idx] else 0
#         for i in duals[d_idx]:
#             if i < m:
#                 intcpt += duals[d_idx][i] * rhs[i]
#             elif duals[d_idx][i] > 1e-5:
#                 intcpt += duals[d_idx][i] * t['lb'][i - m]
#             elif duals[d_idx][i] < -1e-5:
#                 intcpt += duals[d_idx][i] * t['ub'][i - m]
#         print("leaf %d : dual sol %d : (%f, %f)" % (t_idx, d_idx, zeta, intcpt))
        