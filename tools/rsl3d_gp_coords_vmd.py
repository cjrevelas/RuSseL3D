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
            nodeId = int(splitted_line[0])
            #vmdFile.write(nodeId + '\n')


            meshFile = open("in.mesh", 'r')
            ii = 0
            for line_counter_3 in meshFile:
                if ("Mesh point coordinates" in line_counter_3):

                    for line_counter_4 in meshFile:
                        ii += 1
                        if (ii == nodeId):
                            splitted_line_2 = line_counter_4.split()
                            xx_nodeId = float(splitted_line_2[0])
                            yy_nodeId = float(splitted_line_2[1])
                            zz_nodeId = float(splitted_line_2[2])


                            vmdFile.write(str(kk) + space + '0' + space + str(xx_nodeId) + space + str(yy_nodeId) + space + str(zz_nodeId) + space + str(nodeId) + '\n')
            meshFile.close()
    else:
        if ("NUMBER OF ATOMS" in line_counter_1):
            vmdFile.write(line_counter_1)
            vmdFile.write(str(int(gnodesFile.readline())+1)+'\n')
        else:
            vmdFile.write(line_counter_1)

vmdFile.write(str(kk+1) + space + '1' + space + str(0.0) + space + str(0.0) + space + str(0.0) + space + '1\n')
exit()
