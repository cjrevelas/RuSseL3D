try:
    gnodesFile = open("in.gnodes", 'r')
except:
    print("ERROR OPENING GNODES FILE in.gnodes")
    exit()

try:
    vmdFile = open("vmd.lammpstrj", 'w')
except:
    print("ERROR OPENING VMD FILE vmd.lammpstrj")

space = "      "
kk = 0
for line_counter_1 in gnodesFile:
    if ("xu yu zu" in line_counter_1):
        vmdFile.write(line_counter_1.rstrip() + ' ' + "origID\n")
        for line_counter_2 in gnodesFile:
            kk += 1
            splitted_line = line_counter_2.split()
            vmdFile.write(str(kk) + space + '0' + space + splitted_line[3] + space + splitted_line[4] + space + splitted_line[5] + space + splitted_line[0] + '\n')
    else:
        if ("NUMBER OF ATOMS" in line_counter_1):
            vmdFile.write(line_counter_1)
            vmdFile.write(str(int(gnodesFile.readline()))+'\n')
        else:
            vmdFile.write(line_counter_1)

exit()
