#------------------------------------------------------------------------------------------------------------#
import os
import sys
import ast
import numpy  as np
import math   as m
import copy   as cp
import random as rnd
import subprocess
import csv
import fileinput
#------------------------------------------------------------------------------------------------------------#
NULLVAL = -1
N_MAX   = 150
N_AVOG  = 6.022 * 10**23

export_thermo    = True
export_phi_smear = False
restart          = False
with_mpi         = True

read_field = 1
init_iter  = 1
num_iter   = 2

compute_dens_profs_every       = 1
compute_indiv_dens_profs_every = 1
compute_field_every            = 1
compute_binary_field_every     = 0
compute_propagators_every      = 0
compute_brush_every            = 1
compute_chainshape_every       = 1
compute_ads_free_every         = 1
compute_chain_ends_every       = 1

if with_mpi:
    exec_file = "home/cjrevelas/bin/RuSseL3D_2023-01-29_MPI"
else:
    exec_file = "home/cjrevelas/bin/RuSseL3D_2023-01-29_SERIAL"

submit_job = "qsub /home/cjrevelas/bin/runqueue.sh"
input_file = "in.input"
log_file   = "o.log"
#------------------------------------------------------------------------------------------------------------#
def run_qsub(directory):
   os.chdir(directory)
   os.system(submit_job)
   os.chdir("..")

   return
#------------------------------------------------------------------------------------------------------------#
def get_directories():
    tempDirsList = os.listdir(path='.')
    dirs = []
    for directory in tempDirsList:
        path = directory + '/' + log_file
        if os.path.exists(path):
            dirs.append(directory)

    return dirs
#------------------------------------------------------------------------------------------------------------#
def modify_input_parameters(input_path):
    standard_entry_length = 16

    if os.path.exists(input_path):
        input = open(input_path, 'r+')

        while True:
            line = input.readline()
            if not line: break
            if "# init field" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(read_field) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# init iter" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(init_iter) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# num iter" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(num_iter) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export dens profs" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_dens_profs_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export indiv dens profs" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_indiv_dens_profs_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export field" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_field_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export binary field" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_binary_field_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export propagators" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_propagators_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export brush thickness" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_brush_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export chains per area profs" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_chainshape_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export ads vs free profs" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_ads_free_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
            if "# export chain ends profs" in line:
                oldLine      = line.rstrip('\n')
                newLine      = oldLine.replace(oldLine.split()[0] + ' ', str(compute_chain_ends_every) + ' ')
                newLineSplit = newLine.split()
                for kk in range(standard_entry_length-len(newLine.split()[0])):
                    newLineSplit[0] += ' '
                cmd = "sed -i \"s|" + oldLine + '|' + ' '.join(newLineSplit) + "|g\" " + input_path
                os.system(cmd)
        input.close()

    return
#------------------------------------------------------------------------------------------------------------#
def restart_calculations(dirs):
    for directory in dirs:
        path = directory + "/" + input_file
        modify_input_parameters(path)
        run_qsub(directory)

    return
#------------------------------------------------------------------------------------------------------------#
def get_box_dimensions(directory):
    path = directory + '/' + log_file

    box_length = []

    if os.path.exists(path):
        log = open(path, 'r')

        while True:
            line = log.readline()

            if "box_length" in line:
                for ii in range(3):
                    line = log.readline().split()
                    box_length.append(float(line[1]))
                break

    return box_length
#------------------------------------------------------------------------------------------------------------#
def cast_fort_float(val):
    return float(val.replace('D', 'E'))
#------------------------------------------------------------------------------------------------------------#
def get_sim_full_path(directory):
    cmd = "readlink -f " + directory
    output = subprocess.getoutput(cmd)

    return output
#------------------------------------------------------------------------------------------------------------#
def get_brush_thickness(directory, filename):
    [mean, stdev, all] = [NULLVAL]*3
    path = directory + '/' + filename

    if os.path.exists(path):
        brush_thickness_file = open(path, 'r')
        while True:
            line = brush_thickness_file.readline()
            if not line: break
            if "mean"  in line: mean  = line.split()[1]
            if "stdev" in line: stdev = line.split()[1]
            if "all"   in line: all   = line.split()[1]

    return [mean, stdev, all]
#------------------------------------------------------------------------------------------------------------#
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

    #for equidistant grafting
    if kind==NULLVAL:
        kind = "uniform"

    return [kind, str(numPoles), str(param1), str(param2), str(index)]
#------------------------------------------------------------------------------------------------------------#
def get_input_params(directory):
    path = directory + '/' + input_file

    sphere_args = []
    face_args   = []

    wall_dist   = NULLVAL
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
        input = open(path, 'r')
        while True:
            line = input.readline()
            if not line: break
            if "# wall dist" in line: wall_dist = cast_fort_float(line.split()[0])
            if "# num nanop" in line:
                n_spheres = int(line.split()[0])
                for ii in range(n_spheres):
                    line = input.readline().split()
                    sphere_args.append(line)
            if "# num faces"  in line:
                n_faces = int(line.split()[0])
                for ii in range(n_faces):
                    line = input.readline().split()
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

    return [use_mx, use_gr, n_spheres, N_mx, N_gr, wall_dist, eos_type, fraction, ds_ed_mx, ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr, xs_crit_gr, sphere_args, face_args]
#------------------------------------------------------------------------------------------------------------#
def get_energies(directory):
    energies = [NULLVAL]*6

    filename = "o.energy_terms"
    path = directory + '/' + filename

    if os.path.exists(path):
        energies_file  = open(path, 'r')
        line = energies_file.readline()
        line = energies_file.readline()
        energies = line.split()
        energies_file.close()

    energies = [energies[ii] for ii in range(6)]

    return energies
#------------------------------------------------------------------------------------------------------------#
def is_finished(directory):
    run_state = -1

    path = directory + '/' + log_file

    if "SUMMARIZED RESULTS" in open(path).read():
        run_state = 1
    else:
        run_state = 0

    return run_state
#------------------------------------------------------------------------------------------------------------#
def get_last_thermo(directory):
    [step, n_gr_chains, max_error, std_error] = [NULLVAL]*4

    path = directory + '/' + log_file

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
#------------------------------------------------------------------------------------------------------------#
def get_phi_smear(directory):
    filename = "o.phi_smear_np1"
    path     = directory + '/' + filename
    rr       = [0.0]*N_MAX
    phi_mx   = [0.0]*N_MAX
    phi_gr   = [0.0]*N_MAX

    try:
        phi_smeared_file = open(path, 'r')
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
        print("FAILED: " + directory)
        pass

    return [rr, phi_mx, phi_gr]
#------------------------------------------------------------------------------------------------------------#
dirs = get_directories()

if export_thermo:
    csvFile   = open("RuSseL3D.csv", 'w')
    csvWriter = csv.writer(csvFile)

    csvWriter.writerow(["dirname", "run_state", "kind", "numPoles", "param1", "param2", "index",           \
                        "r_np_eff (A)", "use_mx", "use_gr", "n_spheres", "N_mx", "N_gr", "n_gr_chains", "interf_area (A2)",\
                        "gdens (A-2)", "free_energy (mJ/m2)", "free_energy (kJ/mol)", "term1 (mJ/m2)", "term2 (mJ/m2)",    \
                        "term3 (mJ/m2)", "term4 (mJ/m2)", "term4_norm (mJ/m2)", "hh_mean (A)", "hh_std (A)", "hh_all (A)", \
                        "hh99_mean (A)", "hh99_std (A)", "hh99_all (A)", "wall_dist (A)", "eos_type", "fraction", "step",  \
                        "max_error", "std_error", "ds_ed_mx", "ds_conv_mx", "xs_crit_mx", "ds_ed_gr", "ds_conv_gr",        \
                        "xs_crit_gr", "sphere_args", "face_args"])

    for directory in dirs:
        [use_mx, use_gr, n_spheres, N_mx, N_gr, wall_dist, eos_type,      \
         fraction, ds_ed_mx, ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr,\
         xs_crit_gr, sphere_args, face_args] = get_input_params(directory)

        [step, n_gr_chains, max_error, std_error] = get_last_thermo(directory)

        box = get_box_dimensions(directory)

        if sphere_args != []:
            r_np_eff    = float(sphere_args[0][1])
            interf_area = n_spheres * 4.0 * np.pi * (r_np_eff - wall_dist)**2
        else:
            r_np_eff    = -1
            interf_area = box[0] * box[1]

        N_mx  = float(N_mx)
        N_gr  = float(N_gr)
        gdens = float(n_gr_chains) / interf_area

        [term1, term2, term3, term4, term4_norm, free_energy] = get_energies(directory)

        free_energy_kJ_mol = (float(free_energy) * 1e-6) * (interf_area * 1e-20) * N_AVOG

        [hh_mean, hh_std, hh_all]       = get_brush_thickness(directory, "o.brush_np1")
        [hh99_mean, hh99_std, hh99_all] = get_brush_thickness(directory, "o.brush99_np1")

        run_state = is_finished(directory)

        full_path = get_sim_full_path(directory)

        csvWriter.writerow([directory, run_state] + get_params_from_dirname(directory) +                    \
                           [r_np_eff, use_mx, use_gr, n_spheres, N_mx, N_gr, n_gr_chains, interf_area,      \
                            gdens, free_energy, free_energy_kJ_mol, term1, term2, term3, term4, term4_norm, \
                            hh_mean, hh_std, hh_all, hh99_mean, hh99_std, hh99_all, wall_dist,              \
                            eos_type, fraction, step, max_error, std_error, ds_ed_mx,                       \
                            ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr, xs_crit_gr,                       \
                            sphere_args, face_args])

    csvFile.close()
#------------------------------------------------------------------------------------------------------------#
if export_phi_smear:
    prof_smeared = {}
    for directory in dirs:
        [use_mx, use_gr, n_spheres, N_mx, N_gr, wall_dist, eos_type, fraction, \
         ds_ed_mx, ds_conv_mx, xs_crit_mx, ds_ed_gr, ds_conv_gr, xs_crit_gr,   \
         sphere_args, face_args] = get_input_params(directory)

        [step, n_gr_chains, max_error, std_error] = get_last_thermo(directory)

        r_np_eff    = float(sphere_args[0][1])
        N_mx        = float(N_mx)
        N_gr        = float(N_gr)
        interf_area = n_spheres * 4.0 * np.pi * (r_np_eff - wall_dist)**2
        gdens       = float(n_gr_chains) / interf_area

        [rr, prof_mx, prof_gr] = get_phi_smear(directory)

        tag = str(r_np_eff) + '_' + str(N_mx) + '_' + str(N_gr) + '_' + str(gdens)
        prof_smeared[tag] = [rr, prof_mx, prof_gr]

    prof_mx_out = open("o.phi_smear_mx.txt", 'w')
    prof_mx_out.write("bin ")

    for tag in prof_smeared:
        prof_mx_out.write("%s " %(tag))
    prof_mx_out.write('\n')

    for binn in range(len(prof_smeared[tag][0])):
        prof_mx_out.write("%d " %(binn))
        for tag in prof_smeared:
            prof_mx_out.write("%f " %(prof_smeared[tag][1][binn]))
        prof_mx_out.write('\n')
    prof_mx_out.close()

    prof_gr_out = open("o.phi_smear_gr.txt", 'w')
    prof_gr_out.write("bin ")

    for tag in prof_smeared:
        prof_gr_out.write("%s " %(tag))
    prof_gr_out.write('\n')

    for binn in range(len(prof_smeared[tag][0])):
        prof_gr_out.write("%d " % (binn) )
        for tag in prof_smeared:
            prof_gr_out.write("%f " %(prof_smeared[tag][2][binn]))
        prof_gr_out.write('\n')
    prof_gr_out.close()
#------------------------------------------------------------------------------------------------------------#
if restart:
    restart_calculations(dirs)

exit()
#------------------------------------------------------------------------------------------------------------#
