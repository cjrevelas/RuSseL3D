# python script to generate vtu files from rsl3d profiles
import os
import sys
import numpy as np

try:
    mesh_file = open("in.mesh", 'r')
except:
    print("ERROR OPENING MESH FILE in.mesh")
    exit()

try:
    solution_file = open("o.phi", 'r')
except:
    print("ERROR OPENING SOLUTION FILE o.phi")

try:
    vtu_file = open("o.vtu", 'w')
except:
    print("ERROR OPENING VTU FILE o.vtu")

for line in mesh_file:
    if ("sdim" in line):
        ndm = int(line.split()[0])

    if ("number of mesh points" in line):
        numnp = int(line.split()[0])
        meshpoint     = np.zeros((ndm, numnp))
        meshpoint_vtu = np.zeros((ndm,numnp))
        solution      = np.zeros(numnp)

    if ("Mesh point coordinates" in line):
        for jj in range(0,numnp):
            lineBuffer = mesh_file.readline()
            for ii in range(0,ndm):
                meshpoint[ii,jj] = float(lineBuffer.split()[ii])

    if ("4 # number of nodes per element" in line):
        nen = int(line.split()[0])
        numel = int(mesh_file.readline().split()[0])

        global_node_id = np.zeros((nen,numel), dtype=int)

        mesh_file.readline()

        for ii in range(0,numel):
            lineBuffer = mesh_file.readline()

            for jj in range(0,nen):
                global_node_id[jj,ii] = lineBuffer.split()[jj]


for line in solution_file:
    if ("% X" in line):
        for jj in range(0,numnp):
            lineBuffer = solution_file.readline()
            solution[jj] = float(lineBuffer.split()[3])

            for ii in range(0,ndm):
                meshpoint_vtu[ii,jj] = float(lineBuffer.split()[ii])

line_1 = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
line_2 = "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
line_3 = "<Piece NumberOfPoints=\"" + str(numnp) + "\" NumberOfCells=\"" + str(numel) + "\">\n"
line_4 = "<DataArray type=\"Float64\" Name=\"Dependent_variable_u\" Format=\"ascii\">\n"
line_5 = "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n"
line_6 = "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n"
line_7 = "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n"
line_8 = "<DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n"

beginUnstructuredGrid = "<UnstructuredGrid>\n"
endUnstructuredGrid = "</UnstructuredGrid>\n"

beginPointData = "<PointData>\n"
endPointData = "</PointData>\n"

beginCellData = "<CellData/>\n"
endCellData = "</CellData/>\n"

beginPoints = "<Points>\n"
endPoints = "</Points>\n"

beginCells = "<Cells>\n"
endCells = "</Cells>\n"

endDataArray = "</DataArray>\n"
endPiece = "</Piece>\n"
endVTKfile = "</VTKFile>"

vtu_file.write(line_1)
vtu_file.write(line_2)
vtu_file.write(beginUnstructuredGrid)
vtu_file.write(line_3)
vtu_file.write('\n')

vtu_file.write(beginPointData)
vtu_file.write(line_4)
# Insert the solution vector
for ii in range(0,numnp):
    vtu_file.write(str(solution[ii]) + '\n')
vtu_file.write(endDataArray)
vtu_file.write(endPointData)
vtu_file.write('\n')

vtu_file.write(beginCellData)
vtu_file.write(beginPoints)
vtu_file.write(line_5)
# Insert the nodal point coordinates
for jj in range(0,numnp):
    for ii in  range(0,ndm):
        #vtu_file.write("%.10f" %(meshpoint[ii,jj]) + ' ')
        vtu_file.write(str(meshpoint[ii,jj]) + ' ')
        #vtu_file.write(str(meshpoint_vtu[ii,jj]) + ' ')
    vtu_file.write('\n')
vtu_file.write(endDataArray)
vtu_file.write(endPoints)
vtu_file.write('\n')

vtu_file.write(beginCells)
vtu_file.write(line_6)
# Insert element connectivity
for ii in range(0,numel):
    for jj in range(0,nen):
        vtu_file.write(str(global_node_id[jj,ii]) + ' ')
    vtu_file.write('\n')
vtu_file.write(endDataArray)
vtu_file.write('\n')

vtu_file.write(line_7)
# Insert array of offsets
kk = 0
offset = 0
for ii in range(0,numel):
    kk += 1
    offset += nen

    vtu_file.write(str(offset) + ' ')
    if (kk%nen==0):
        vtu_file.write('\n')
vtu_file.write('\n')
vtu_file.write(endDataArray)
vtu_file.write('\n')

vtu_file.write(line_8)
# Insert array of types
kk = 0
tet_type = 10
for ii in range(0,numel):
    kk += 1

    vtu_file.write(str(tet_type) + ' ')
    if (kk%(nen*nen)==0):
        vtu_file.write('\n')
vtu_file.write('\n')
vtu_file.write(endDataArray)
vtu_file.write(endCells)
vtu_file.write('\n')

vtu_file.write(endPiece)
vtu_file.write(endUnstructuredGrid)
vtu_file.write(endVTKfile)

mesh_file.close()
solution_file.close()
vtu_file.close()

exit()
