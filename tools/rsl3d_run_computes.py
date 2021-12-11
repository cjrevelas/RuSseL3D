import os
import sys
import numpy  as np
import math   as m
import copy   as cp
import subprocess
import fileinput

with_mpi = True

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
    exec_file = "/home/cjrevelas/bin/RuSseL3D_2021-01-06_MPI"
else:
    exec_file = "/home/cjrevelas/bin/RuSseL3D_2021-01-06_SERIAL"

submit_job = "qsub /home/cjrevelas/bin/runqueue.sh" # + exec_file 
input_name = "input.in.txt"
log_name   = "o.log.out.txt"


def run_qsub(dir):
    os.chdir(dir)
    os.system(submit_job)
    os.chdir("..")

    return


def modify_input_params(input_path):
    standard_entry_length = 16

    if os.path.exists(input_path):
        input_file = open(input_path, 'r+')

        while True:
            line = input_file.readline()
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
        input_file.close() 

    return 


def restart_calculations():
    tempDirsList = os.listdir(path='.')
    dirs = []
    for dir in tempDirsList:
        path = dir + '/' + log_name
        if os.path.exists(path):
            dirs.append(dir)

    for dir in dirs:
        path = dir + '/' + input_name
        modify_input_params(path)
        run_qsub(dir)

    return


restart_calculations()
exit()
