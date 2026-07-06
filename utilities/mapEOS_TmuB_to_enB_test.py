#!/usr/bin/env python3
# Copyright Chun Shen @ 2018
# Updated to replace deprecated scipy.interpolate.interp2d with
# scipy.interpolate.RegularGridInterpolator (SciPy >= 1.14)

from numpy import *
import sys
from os import path
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import multiprocessing as mp

g1 = 2.*(3.*3.-1) + 7.*3.*2.5/2.
pi = 3.14159265358979323846264338327950288419716939937510

try:
    EOS_table_folder = str(sys.argv[1])
except:
    print("Usage: python3 mapEOS_TmuB_to_enB.py EOS_table_folder")
    exit(1)

ACCURACY = 1e-6
MAXITER  = 100
hbarC    = 0.19733

filename_string = EOS_table_folder
#filename_string = EOS_table_folder.split("Files_")[1].split("/")[0].split("_SN")[0]

print(filename_string)

ed_table = loadtxt(path.join(EOS_table_folder,
                             "EnerDens.dat"))    # e/T^4
nB_table = loadtxt(path.join(EOS_table_folder,
                             "BarDens.dat"))     # nB/T^3
P_table  = loadtxt(path.join(EOS_table_folder,
                             "Press.dat"))       # P/T^4
CorrLength_table  = loadtxt(path.join(EOS_table_folder,
                             "CorrLength.dat"))  # fm
Cs2_table  = loadtxt(path.join(EOS_table_folder,
                             "Cs2.dat"))
Chi2_table  = loadtxt(path.join(EOS_table_folder,
                             "Chi2.dat"))        # chi2/T^2
s_table  = loadtxt(path.join(EOS_table_folder,
                             "EntrDens.dat"))    # s/T^3


#n_muB = 471
#muB = linspace(0.000, 0.470, n_muB)
n_muB = 601
muB = linspace(0.000, 0.600, n_muB)
#n_muB = 451
#muB = linspace(0.000, 0.450, n_muB)
n_T   = 771
T   = linspace(0.030, 0.800, n_T)

muB_table, T_table = meshgrid(muB, T)
ed_table = ed_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
nB_table = nB_table[:, 2].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)        # 1/fm^3
P_table  = P_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)         # GeV/fm^3
Cs2_table  = Cs2_table[:, 2].reshape(n_T, n_muB)
CorrLength_table = CorrLength_table[:, 2].reshape(n_T, n_muB)
Chi2_table = Chi2_table[:, 2].reshape(n_T, n_muB)*(T_table**2.)/(hbarC**3.)    # 1/(GeV*fm^3)
s_table  = s_table[:, 2].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)  # 1/fm^3



print("e_min = %.5e GeV/fm^3, e_max = %.5e GeV/fm^3"
        % (matrix(ed_table).min(), matrix(ed_table).max()))
print("nB_min = %.5e 1/fm^3, nB_max = %.5e 1/fm^3"
        % (matrix(nB_table).min(), matrix(nB_table).max()))

# NOTE:
# RegularGridInterpolator expects the `values` array shaped in the
# same order as the grid tuple. We want the grid to be (muB, T) so
# `values` must have shape (len(muB), len(T)). Current tables are
# shaped (n_T, n_muB) -> (len(T), len(muB)) so we take the transpose.

p_interp    = RegularGridInterpolator((muB, T), P_table.T,  method='linear', bounds_error=False, fill_value=None)
e_interp    = RegularGridInterpolator((muB, T), ed_table.T, method='linear', bounds_error=False, fill_value=None)
nB_interp   = RegularGridInterpolator((muB, T), nB_table.T, method='linear', bounds_error=False, fill_value=None)
cs2_interp  = RegularGridInterpolator((muB, T), Cs2_table.T, method='linear', bounds_error=False, fill_value=None)
corr_interp = RegularGridInterpolator((muB, T), CorrLength_table.T, method='linear', bounds_error=False, fill_value=None)
chi2_interp = RegularGridInterpolator((muB, T), Chi2_table.T, method='linear', bounds_error=False, fill_value=None)
s_interp    = RegularGridInterpolator((muB, T), s_table.T, method='linear', bounds_error=False, fill_value=None)

# wrapper functions to mimic old interp2d call signature f(x,y)
def f_p(muB_val, T_val):
    return p_interp([[muB_val, T_val]])[0]

def f_e(muB_val, T_val):
    return e_interp([[muB_val, T_val]])[0]

def f_nB(muB_val, T_val):
    return nB_interp([[muB_val, T_val]])[0]

def f_cs2(muB_val, T_val):
    return cs2_interp([[muB_val, T_val]])[0]

def f_corr(muB_val, T_val):
    return corr_interp([[muB_val, T_val]])[0]

def f_chi2(muB_val, T_val):
    return chi2_interp([[muB_val, T_val]])[0]

def f_s(muB_val, T_val):
    return s_interp([[muB_val, T_val]])[0]


def smoothstep01(x):
    x = min(max(x, 0.0), 1.0)
    return x*x*x*(10.0 - 15.0*x + 6.0*x*x)

def safe_s(e, p, muB, nB, T):
    s_val = (e + p - muB*nB)/max(T, 1e-12)
    if not isfinite(s_val) or s_val < 0.0:
        return 0.0
    return s_val

def eval_at(muB_val, T_val, ed_local=None, nB_local=None):
    P_local     = f_p(muB_val, T_val)
    cs2_local   = f_cs2(muB_val, T_val)
    corrL_local = f_corr(muB_val, T_val)
    chi2_local  = f_chi2(muB_val, T_val)

    if ed_local is not None and nB_local is not None:
        s_local = safe_s(ed_local, P_local, muB_val, nB_local, T_val)
    else:
        s_local = f_s(muB_val, T_val)

    return P_local, T_val, muB_val, cs2_local, corrL_local, chi2_local, s_local

def binary_search_1d(ed_local, muB_local):
    T_low = 0.030
    T_high = 0.800

    e_low = f_e(muB_local, T_low)
    e_high = f_e(muB_local, T_high)

    if not isfinite(e_low) or not isfinite(e_high):
        return T_low, False

    if ed_local <= e_low:
        return T_low, False

    if ed_local >= e_high:
        return T_high, False

    for _ in range(MAXITER):
        T_mid = 0.5*(T_low + T_high)
        e_mid = f_e(muB_local, T_mid)

        if not isfinite(e_mid):
            return T_mid, False

        if ed_local < e_mid:
            T_high = T_mid
        else:
            T_low = T_mid

        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)

        if rel_err < ACCURACY or abs_err < ACCURACY*1e-2:
            return T_mid, True

    return T_mid, True

NB_TOL = 1e-12

def invert_inside_image(ed_local, nB_local):
    muB_low = 0.0
    muB_high = muB[-1]

    # Caso especial: nB = 0 deve cair em muB = 0
    if abs(nB_local) < NB_TOL:
        T0, valid_T0 = binary_search_1d(ed_local, 0.0)

        if valid_T0:
            return T0, 0.0, True, 0.0, f_nB(muB_high, binary_search_1d(ed_local, muB_high)[0])

        # se a energia nem existe na tabela em muB=0, aí sim é inválido
        return T0, 0.0, False, 0.0, 0.0

    T_low_mu, valid_low = binary_search_1d(ed_local, muB_low)
    T_high_mu, valid_high = binary_search_1d(ed_local, muB_high)

    nB_min_phys = f_nB(muB_low, T_low_mu)
    nB_max_phys = f_nB(muB_high, T_high_mu)

    # força simetria em muB=0
    if abs(nB_min_phys) < 1e-10:
        nB_min_phys = 0.0

    if nB_local < nB_min_phys:
        return T_low_mu, muB_low, False, nB_min_phys, nB_max_phys

    if nB_local > nB_max_phys:
        return T_high_mu, muB_high, False, nB_min_phys, nB_max_phys

    lo = muB_low
    hi = muB_high

    for _ in range(MAXITER):
        mid = 0.5*(lo + hi)
        T_mid, valid_T = binary_search_1d(ed_local, mid)
        nB_mid = f_nB(mid, T_mid)

        if not isfinite(nB_mid):
            return T_mid, mid, False, nB_min_phys, nB_max_phys

        if nB_local < nB_mid:
            hi = mid
        else:
            lo = mid

        abs_err = abs(nB_mid - nB_local)
        rel_err = abs_err/abs(nB_mid + nB_local + 1e-15)

        if rel_err < ACCURACY or abs_err < ACCURACY*1e-2:
            return T_mid, mid, True, nB_min_phys, nB_max_phys

    return T_mid, mid, True, nB_min_phys, nB_max_phys

def eos_conformal(ed_local):
    e_safe = max(ed_local, 1e-14)

    p_conf   = e_safe/3.0
    t_conf   = (30.0*e_safe/(pi*pi*g1))**0.25
    s_conf   = (30.0*e_safe/(pi*pi*g1))**0.75 * 2.0*pi*pi*g1/45.0
    cs2_conf = 1.0/3.0

    return p_conf, t_conf, 0.0, cs2_conf, 0.0, 0.0, s_conf


def complete_outside_image(ed_local, nB_local,
                           T_edge, muB_edge,
                           nB_min_phys, nB_max_phys,
                           drhob):

    P_local     = f_p(muB_edge, T_edge)
    cs2_local   = f_cs2(muB_edge, T_edge)
    corr_local  = f_corr(muB_edge, T_edge)
    chi2_local  = f_chi2(muB_edge, T_edge)

    # mantém T e muB na borda física
    T_local   = T_edge
    muB_local = muB_edge

    # entropia consistente
    s_local = (ed_local + P_local - muB_local*nB_local)/max(T_local, 1e-12)

    if not isfinite(s_local) or s_local < 0.0:
        s_local = f_s(muB_edge, T_edge)

    return P_local, T_local, muB_local, cs2_local, corr_local, chi2_local, s_local


def invert_or_complete(ed_local, nB_local, drhob):

    T_local, muB_local, valid, nB_min_phys, nB_max_phys = \
        invert_inside_image(ed_local, nB_local)

    if valid:
        P_local, T_local, muB_local, cs2_local, corrL_local, chi2_local, s_local = \
            eval_at(muB_local, T_local, ed_local, nB_local)

        return (P_local, T_local, muB_local,
                cs2_local, corrL_local, chi2_local,
                s_local, 1.0)

    if nB_local > nB_max_phys:
        T_edge, _ = binary_search_1d(ed_local, muB[-1])
        muB_edge = muB[-1]
    else:
        T_edge, _ = binary_search_1d(ed_local, 0.0)
        muB_edge = 0.0

    P_local, T_local, muB_local, cs2_local, corrL_local, chi2_local, s_local = \
        complete_outside_image(
            ed_local, nB_local,
            T_edge, muB_edge,
            nB_min_phys, nB_max_phys,
            drhob
        )

    return (P_local, T_local, muB_local,
            cs2_local, corrL_local, chi2_local,
            s_local, 0.0)

#T_local, muB_local = binary_search_2d(1.0, 0.02)
#print(T_local, muB_local, f_e(muB_local, T_local), f_nB(muB_local, T_local))


#ed_bounds = [0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 219.355, 2209.36]
#ne_list = [14, 21, 32, 43, 201, 201, 201]
#nB_bounds = [0.00499, 0.01495, 0.04475, 0.498, 3.49, 12.45, 39.8]
#nnB_list  = [501, 301, 181, 251, 351, 251, 201]


ed_bounds = [0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 80, 650.]
ne_list   = [61, 60, 61, 122, 200, 400, 400]
#nB_bounds = [0.0025, 0.015, 0.045, 0.25, 3.5, 9, 14]
nB_bounds = [0.0025, 0.015, 0.045, 0.25, 3.5, 8.0, 12]
nnB_list  = [501, 301, 181, 251, 351, 251, 201]
'''
ed_bounds = [0.0, 0.02, 0.12, 0.43, 1.5, 6.0, 20.0, 60.0]
ne_list   = [100, 120, 140, 180, 220, 220, 140]
nB_bounds = [0.01, 0.05, 0.3, 1.5, 5.0, 10.0, 17.5]
nnB_list  = [401, 301, 251, 251, 221, 201, 161]

'''

'''
ed_bounds = [0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 650.]
ne_list   = [61, 60, 61, 122, 200, 400]
nB_bounds = [0.0025, 0.015, 0.045, 0.25, 3.5, 12.0]
nnB_list  = [501, 301, 181, 251, 351, 251]
'''



# generate tables
for itable in range(len(ne_list)):
    print("Generating table {0:d} ({1:d} x {2:d})... ".format(
        itable, ne_list[itable], nnB_list[itable]))

    # generate the grid arrays
    ed_list  = linspace(ed_bounds[itable], ed_bounds[itable+1], ne_list[itable])
    de       = ed_list[1] - ed_list[0]
    nB_list  = linspace(0.0, nB_bounds[itable], nnB_list[itable])
    drhob    = nB_list[1] - nB_list[0]
    p_list   = zeros([len(ed_list), len(nB_list)])
    T_list   = zeros([len(ed_list), len(nB_list)])
    muB_list = zeros([len(ed_list), len(nB_list)])
    cs2_list = zeros([len(ed_list), len(nB_list)])
    corr_list = zeros([len(ed_list), len(nB_list)])
    chi2_list = zeros([len(ed_list), len(nB_list)])
    s_list   = zeros([len(ed_list), len(nB_list)])
    valid_list = zeros([len(ed_list), len(nB_list)])

    def invert_EOS_tables(idx):
        ie       = int(idx/len(nB_list))
        inB      = idx % len(nB_list)
        ed_local = ed_list[ie]
        nB_local = nB_list[inB]

        P_local, T_local, muB_local, cs2_local, corrL_local, chi2_local, s_local, valid_flag = \
            invert_or_complete(ed_local, nB_local, drhob)

        return (ie, inB, P_local, T_local, muB_local,
                cs2_local, corrL_local, chi2_local, s_local, valid_flag)

    # calculate every points in parallel
    inputs  = range(len(ed_list)*len(nB_list))
    pool    = mp.Pool(processes=70)
    results = pool.map(invert_EOS_tables, inputs)
    pool.close()
    pool.join()
    #results = []
    #for idx_ in inputs:
    #    results.append(invert_EOS_tables(idx_))

    # fill into the tables
    for idx in range(len(results)):
        ie                = results[idx][0]
        inB               = results[idx][1]
        p_list[ie, inB]   = results[idx][2]
        T_list[ie, inB]   = results[idx][3]
        muB_list[ie, inB] = results[idx][4]
        cs2_list[ie, inB] = results[idx][5]
        corr_list[ie, inB] = results[idx][6]
        chi2_list[ie, inB] = results[idx][7]
        s_list[ie, inB]    = results[idx][8]
        valid_list[ie, inB]    = results[idx][9]

    # save to files
    savetxt(EOS_table_folder+"/BEST_eos_p_{0:d}.dat".format(itable), p_list,
            fmt='%.8e', delimiter="  ",
            header="%.8e %.8e %d %.8e %.8e %d"
                    % (ed_bounds[itable], de, ne_list[itable],
                    0.0, drhob, nnB_list[itable]),
            comments='')
    savetxt(EOS_table_folder+"/BEST_eos_T_{}.dat".format(itable), T_list,
            fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_muB_{}.dat".format(itable), muB_list,
            fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_cs2_{}.dat".format(itable), cs2_list,
            fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_corrLength_{}.dat".format(itable), corr_list,
            fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_chi2_{}.dat".format(itable), chi2_list,
            fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_s_{}.dat".format(itable), s_list,
        fmt='%.8e', delimiter="  ")
    savetxt(EOS_table_folder+"/BEST_eos_valid_{}.dat".format(itable),
        valid_list, fmt='%.1f', delimiter="  ")