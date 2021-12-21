import os
import sys
import ast
import numpy as np
import math as m
import copy as cp
import random as rnd

input_template = "in.input.template3d"
input_name     = "in.input"

if with_mpi:
    exec_file = "/home/cjrevelas/bin/RuSseL3D_2021-01-06_MPI"
else:
    exec_file = "/home/cjrevelas/bin/RuSseL3D_2021-01-06_SERIAL"

submit_job     = "qsub /home/cjrevelas/bin/runqueue.sh" + exec_file


def get_fraction_from_chainlen(chain_length):

    return 0.17 / float(chain_length)           # the numerator here may needs to be updated for RuSseL3D


def refresh_input():
    os.system("cp " + input_template + ' ' + input_name)

    return


def set_int_flag(flag_type, value):
    strval = str(value)
    for ii in range(len(flag_type)-len(strval)):
        strval += ' '
    cmd = "sed -i 's/" + flag_type + '/' + strval + "/g' " + input_name
    os.system(cmd)

    return


def set_txt_flag(flag_type, value):
    strval = str(value)
    for ii in range(len(flag_type)-len(strval)):
        strval += ' '
    cmd = "sed -i 's/" + flag_type + '/' + strval + "/g' " + input_name
    os.system(cmd)

    return


def set_float_flag(flag_type, value):
    strval = "{:.9e}".format(value)
    for ii in range(len(flag_type)-len(strval)):
        strval += ' '
    cmd = "sed -i 's/" + flag_type + '/' + strval + "/g' " + input_name
    os.system(cmd)

    return


def run_qsub(path):
    os.chdir(path)
    os.system("qsub " + qsub_path + ' ' + exec_file)
    os.chdir("..")

    return


def generate_simulation(r_np, temp, frac, use_mx, N_mx, use_gr, N_gr, gdens):
    folder_name = "Rnp" + str(round(r_np,0)) + '_' + "frac" + str(round(frac,6)) + '_'

    if use_mx:
        folder_name += "Nmx" + str(N_mx) + '_'

    if use_gr:
        folder_name += "Ngr" + str(N_gr) + '_' + "gdens" + str(round(gdens,5))

    if use_mx:
        use_mx = 1
    else:
        use_mx = 0

    if use_gr:
        use_gr = 1
    else:
        use_gr = 0

    os.system("mkdir " + folder_name)

    refresh_input()

    set_float_flag("FLAG_RADIUS"  , r_np)
    set_float_flag("FLAG_TEMP"    , temp)
    set_float_flag("FLAG_FRACTION", frac)
    set_float_flag("FLAG_GDENS"   , gdens)

    set_txt_flag("FLAG_MATRIX" , use_mx)
    set_txt_flag("FLAG_GRAFTED", use_gr)

    set_int_flag("FLAG_NM", N_mx)
    set_int_flag("FLAG_NG", N_gr)

    #move input in the appropriate folder
    os.system("mv " + input_name + ' ' + folder_name)
    os.system("cp " + mesh       + ' ' + folder_name)
    os.system("cp " + gnodes     + ' ' + folder_name)

    run_qsub(folder_name)

#default flags
temp     = 500.0
r_np     = 80.0
frac = 0.032
N_gr     = 20.0
N_mx     = 40.0
gdens    = 0.0
use_mx   = 1
use_gr   = 1

#generate_simulation(r_np, temp, frac, use_mx, N_mx, use_gr, N_gr, gdens)

'''
for N_mx in [768.0]:
  N_gr = N_mx
  for gdens in [0.001, 0.004, 0.008, 0.012, 0.016]:

     N_max = max(N_mx, N_gr)
     for r_np in [640.0, 1280.0, 2560.0]:
         if gdens == 0.001: lx = 500.0
         if gdens == 0.004: lx = 600.0
         if gdens == 0.008: lx = 750.0
         if gdens == 0.012: lx = 1050.0
         if gdens == 0.016: lx = 1350.0

         frac = get_fraction_from_chainlen(N_max)

         use_mx = 1
         use_gr = 1

         generate_simulation(r_np, temp, frac, use_mx, N_mx, use_gr, N_gr, gdens)
'''
