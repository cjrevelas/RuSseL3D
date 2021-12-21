try:
    gnodesFile = open("in.gnodes", 'r')
except:
    print("ERROR OPENING GNODES FILE in.gnodes")
    exit()

try:
    vmdFile = open("vmd.lammpstrj", 'w')
except:
    print("ERROR OPENING VMD FILE vmd.lammpstrj")

number_of_particles = 2

space = "      "
kk = 0
for line_counter_1 in gnodesFile:
    if ("xu yu zu" in line_counter_1):
        vmdFile.write(line_counter_1.rstrip() + ' ' + "origID\n")
        for line_counter_2 in gnodesFile:
            kk += 1
            splitted_line = line_counter_2.split()
            vmdFile.write(str(kk) + space + '0' + space + splitted_line[2] + space + splitted_line[3] + space + splitted_line[4] + space + splitted_line[0] + '\n')
    else:
        if ("NUMBER OF ATOMS" in line_counter_1):
            vmdFile.write(line_counter_1)
            vmdFile.write(str(int(gnodesFile.readline())+number_of_particles)+'\n')
        else:
            vmdFile.write(line_counter_1)

extra_atom = str(kk+1) + space + '1' + space + "-31.0" + space + "0.0" + space + "0.0" + space + '0' + '\n'
vmdFile.write(extra_atom)
extra_atom = str(kk+2) + space + '1' + space + "+31.0" + space + "0.0" + space + "0.0" + space + '0'
vmdFile.write(extra_atom)

exit()
