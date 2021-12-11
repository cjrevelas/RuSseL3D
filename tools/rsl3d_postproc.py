import os
import sys
import ast
import numpy  as np
import math   as m
import copy   as cp
import random as rnd
import subprocess
import csv 

NULLVAL = -1
N_MAX   = 150
N_AVOG  = 

export_thermo    = True
export_phi_smear = False


def cast_fort_float(val):
    return float(val.replace('D', 'E'))


def get_sim_full_path(dir):
    cmd = "readlink -f " + dir
    output = subprocess.getoutput(cmd)

    return output


def get_brush_thickness(dir, filename):
    [mean, stdev, all] = [NULLVAL]*3
    path = dir + '/' + filename

    if os.path.exists(path):
        brush_thickness_file = open(dir + '/' + filename, 'r')
        while True: 
            line = brush_thickness_file.readline() 
            if not line: break
            if "mean"  in line: mean  = line.split()[1]
            if "stdev" in line: stdev = line.split()[1]
            if "all"   in line: all   = line.split()[1]

    return [mean, stdev, all]


def get_params_from_dirname(dirname):
    dirname_new = dirname.replace('_',' ')
    dirname_splitted = dirname_new.split()

    param1   = 0.0
    param2   = 0
    kind     = NULLVAL
    index    = 1
    numPoles = 0

    if "1pole" in dirname:
        numPoles = 1
    elif "2poles" in dirname:
        numPoles = 2
   
    if "plus" in dirname:
        kind = "attractive"
    elif "minus" in dirname:
        kind = "repulsive"
    elif "homo" in dirname:
        kind = "homo"

    temp1 = ''
    for subString in dirname_splitted:
        if ("plus" in subString) or ("minus" in subString):
            for char in subString:
                if (char.isdigit()) or (char=='.'):
                    temp1 += char
            param1 = float(temp1)
        if subString.isnumeric():
            param2 = int(subString)
        if "index" in subString:
            index = int(subString.replace("index",''))

    #for the uniform grafting point distributions
    if kind==NULLVAL:
        kind = "uniform"

    return [kind, str(numPoles), str(param1), str(param2), str(index)]


def get_input_params(dir):
    filename = "input.in.txt"
    path = dir + '/' + filename

    sphere_args = []
    face_args   = []

    wall_dist   = NULLVAL
    interf_area = NULLVAL
    n_spheres   = NULLVAL
    eos_type    = NULLVAL
    fraction    = NULLVAL
    use_mx      = NULLVAL
    N_mx        = NULLVAL
    ds_ed_mx    = NULLVAL
    ds_conv_mx  = NULLVAL
    xs_crit_mx  = NULLVAL
    use_gr      = NULLVAL
    N_gr        = NULLVAL
    ds_ed_gr    = NULLVAL
    ds_conv_gr  = NULLVAL
    xs_crit_gr  = NULLVAL

    if os.path.exists(path):
        input_file = open(dir + '/' + filename, 'r')
        while True: 
            line = input_file.readline() 
            if not line: break
            if "# wall dist" in line: wall_dist   = cast_fort_float(line.split()[0])
            if "# area"      in line: interf_area = cast_fort_float(line.split()[0])
            if "# num nanop" in line:
                n_spheres = int(line.split()[0])
                for ii in range(n_spheres):
                    line = input_file.readline().split()
                    sphere_args.append(line)
            if "# num faces"  in line:
                n_faces = int(line.split()[0])
                for ii in range(n_faces):
                    line = input_file.readline().split()
                    face_args.append(line)
            if "# eos type"            in line: eos_type = line.split()[0]
            if "# fraction"            in line: fraction = cast_fort_float(line.split()[0])
            if "# use matrix"          in line: use_mx   = line.split()[0]
            if "# chain length matrix" in line: N_mx     = cast_fort_float(line.split()[0])
            if "# contour step matrix" in line:
                ds_ed_mx   = line.split()[0]
                ds_conv_mx = line.split()[1]
            if "# crit contour matrix" in line:
                xs_crit_mx = line.split()[0]
            if "# use grafted"          in line: use_gr = line.split()[0]
            if "# chain length grafted" in line: N_gr   = cast_fort_float(line.split()[0])
            if "# contour step grafted" in line:
                ds_ed_gr   = line.split()[0]
                ds_conv_gr = line.split()[1]
            if "# crit contour grafted"  in line:
                xs_crit_gr = line.split()[0]

    return [use_mx, use_gr, N_mx, N_gr, interf_area, wall_dist, eos_type, fraction, ds_ed_mx, ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr, xs_crit_gr, sphere_args, face_args]


def get_energies(dir):
    energies = [NULLVAL]*6

    filename = "o.energy_terms.out.txt"
    path = dir + '/' + filename

    if os.path.exists(path):
        energies_file  = open(dir + '/' + filename, 'r')
        line = energies_file.readline()
        line = energies_file.readline()
        energies = line.split()
        energies_file.close()

    energies = [energies[ii] for ii in range(6)]

    return energies


def is_finished(dir):
    run_state = -1

    filename = "o.log.out.txt"
    path = dir + '/' + filename

    if "SUMMARIZED RESULTS" in open(path).read():
        run_state = 1
    else:
        run_state = 0

    return run_state


def get_last_thermo(dir):
    [step, n_gr_chains, max_error, std_error] = [NULLVAL]*4

    filename = "o.log.out.txt"

    path = dir + '/' + filename

    if os.path.exists(path):
        for line in reversed(open(path).readlines()):
            linesplit = line.split()
            if linesplit:
                if linesplit[0].isnumeric():
                    step        = linesplit[0]
                    n_gr_chains = linesplit[3]
                    max_error   = linesplit[4]
                    std_error   = linesplit[5]
                    return [step, n_gr_chains, max_error, std_error]

    return [step, n_gr_chains, max_error, std_error]


def get_phi_smear(dir):
    filename = "o.phi_smear_np1.out.txt"
    path     = dir + '/' + filename
    rr       = [0.0]*N_MAX
    phi_mx   = [0.0]*N_MAX
    phi_gr   = [0.0]*N_MAX

    try:
        phi_smeared_file = open(dir + '/' + filename, 'r')
        phi_smeared_file.readline() 
        for ii in range(N_MAX):
            line = phi_smeared_file.readline()
            if not line: break
            line = line.split()
            rr[ii]     = float(line[0])
            phi_mx[ii] = float(line[1])
            phi_gr[ii] = float(line[2])
        phi_smeared_file.close()
    except:
        print("FAILED: " + dir)
        pass

    return [rr, phi_mx, phi_gr]


tempDirsList = os.listdir(path='.')
dirs = []
for dir in tempDirsList:
    filename = "o.log.out.txt"
    path = dir + '/' + filename
    if os.path.exists(path):
        dirs.append(dir)

if export_thermo:
    csvFile   = open("RuSseL3D.csv", 'w')
    csvWriter = csv.writer(csvFile)

    csvWriter.writerow(["dirname", "run_state", "kind", "numPoles", "param1", "param2", "index",           \
                        "r_np_eff", "use_mx", "use_gr", "N_mx", "N_gr", "n_gr_chains", "interf_area",      \
                        "gdens", "free_energy", "term1", "term2", "term3", "term4", "term4_norm",          \
                        "hh_mean", "hh_std", "hh_all", "hh99_mean", "hh99_std", "hh99_all", "wall_dist",   \
                        "eos_type", "fraction", "step", "max_error", "std_error", "ds_ed_mx", "ds_conv_mx",\
                        "xs_crit_mx", "ds_ed_gr", "ds_conv_gr", "xs_crit_gr", "sphere_args", "face_args"])

    for dir in dirs:
        [use_mx, use_gr, N_mx, N_gr, interf_area, wall_dist, eos_type,    \
         fraction, ds_ed_mx, ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr,\
         xs_crit_gr, sphere_args, face_args] = get_input_params(dir)

        [step, n_gr_chains, max_error, std_error] = get_last_thermo(dir)

        r_np_eff = float(sphere_args[0][1])

        N_mx = float(N_mx)
        N_gr = float(N_gr)

        gdens = float(n_gr_chains) / float(interf_area)

        [term1, term2, term3, term4, term4_norm, free_energy] = get_energies(dir)

        [hh_mean, hh_std, hh_all]       = get_brush_thickness(dir, "o.brush_np1.out.txt")
        [hh99_mean, hh99_std, hh99_all] = get_brush_thickness(dir, "o.brush99_np1.out.txt")

        run_state = is_finished(dir)

        full_path = get_sim_full_path(dir)

        csvWriter.writerow([dir, run_state] + get_params_from_dirname(dir) +                  \
                           [r_np_eff, use_mx, use_gr, N_mx, N_gr, n_gr_chains, interf_area,   \
                            gdens, free_energy, term1, term2, term3, term4, term4_norm,       \
                            hh_mean, hh_std, hh_all, hh99_mean, hh99_std, hh99_all, wall_dist,\
                            eos_type, fraction, step, max_error, std_error, ds_ed_mx,         \
                            ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr, xs_crit_gr,         \
                            sphere_args, face_args])
    csvFile.close()
                                        
if export_phi_smear:
    prof_smeared = {}
    for dir in dirs:
        [wall_dist, interf_area, eos_type, fraction, use_mx, N_mx, \
         ds_ed_mx, ds_conv_mx, xs_crit_mx, use_gr, N_gr, ds_ed_gr, \
         ds_conv_gr, xs_crit_gr, sphere_args, face_args] = get_input_params(dir)

        [step, n_gr_chains, max_error, std_error] = get_last_thermo(dir)

        R_np  = float(sphere_args[0][1])
        N_mx  = float(N_mx)
        N_gr  = float(N_gr)
        gdens = float(n_gr_chains) / float(interf_area)

        [rr, prof_mx, prof_gr] = get_phi_smear(dir)

        tag = str(R_np) + '_' + str(N_mx) + '_' + str(N_gr) + '_' + str(gdens)
        prof_smeared[tag] = [rr, prof_mx, prof_gr]

    prof_mx_out = open("o.phi_smear_mx.txt", 'w')
    prof_mx_out.write("bin ")

    for tag in prof_smeared:
        prof_mx_out.write("%s " %(tag))
    prof_mx_out.write('\n')

    for bin in range(len(prof_smeared[tag][0])):
        prof_mx_out.write("%d " %(bin))
        for tag in prof_smeared:
            prof_mx_out.write("%f " %(prof_smeared[tag][1][bin]))
        prof_mx_out.write('\n')
    prof_mx_out.close()

    prof_gr_out = open("o.phi_smear_gr.txt", 'w')
    prof_gr_out.write("bin ")

    for tag in prof_smeared:
        prof_gr_out.write("%s " %(tag))
    prof_gr_out.write('\n')

    for bin in range(len(prof_smeared[tag][0])):
        prof_gr_out.write("%d " % (bin) )
        for tag in prof_smeared:
            prof_gr_out.write("%f " %(prof_smeared[tag][2][bin]))
        prof_gr_out.write('\n')
    prof_gr_out.close()
