import numpy as np
import os
import sys

try:
    mesh_file = open("in.mesh", 'r')
except:
    print("ERROR OPENING INPUT FILE in.mesh")
    exit()

try:
    meshpoints_file = open("o.meshpoints", 'w')
except:
    print("ERROR OPENING OUTPUT FILE o.meshpoints")
    exit()

try:
    elem_con_file_type_zero = open("o.vertex_elements", 'w')
except:
    print("ERROR OPENING OUTPUT FILE o.vertex_elements")
    exit()

try:
    elem_con_file_type_one = open("o.edge_elements", 'w')
except:
    print("ERROR OPENING OUTPUT FILE o.edge_elements")
    exit()

try:
    elem_con_file_type_two = open("o.face_elements", 'w')
except:
    print("ERROR OPENING OUTPUT FILE o.face_elements")
    exit()

try:
    elem_con_file_type_three = open("o.domain_elements", 'w')
except:
    print("ERROR OPENING OUTPUT FILE o.domain_elements")
    exit()

try:
    xface1 = open("o.xface1",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.xface1")
    exit()

try:
    xface2 = open("o.xface2",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.xface2")
    exit()

try:
    yface1 = open("o.yface1",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.yface1")
    exit()

try:
    yface2 = open("o.yface2",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.yface2")
    exit()

try:
    zface1 = open("o.zface1",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.zface1")
    exit()

try:
    zface2 = open("o.zface2",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.zface2")
    exit()

try:
    node_pairing_xx = open("o.node_pairing_xx",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.node_pairing_xx")
    exit()

try:
    node_pairing_yy = open("o.node_pairing_yy",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.node_pairing_yy")
    exit()

try:
    node_pairing_zz = open("o.node_pairing_zz",'w')
except:
    print("ERROR OPENING OUTPUT FILE o.node_pairing_zz")
    exit()

node_hash          = {}
vertex_entity_hash = {}
edge_entity_hash   = {}
face_entity_hash   = {}
domain_entity_hash = {}

tol = 1e-12

for line in mesh_file:
    if ("sdim" in line):
        ndm = int(line.split()[0])

    if ("number of mesh points" in line):
        numnp = int(line.split()[0])
        meshpoint = np.zeros((ndm,numnp))

    if ("Mesh point coordinates" in line):
        for jj in range(0,numnp):
            coords = mesh_file.readline()
            meshpoints_file.write(coords)
            coords = coords.split()
            meshpoint[0,jj] = float(coords[0])
            meshpoint[1,jj] = float(coords[1])
            meshpoint[2,jj] = float(coords[2])

            node_hash[jj+1] = [meshpoint[0,jj], meshpoint[1,jj], meshpoint[2,jj]]

    if ("1 # number of nodes per element" in line):
        nen_type_zero   = int(line.split()[0])
        numel_type_zero = int(mesh_file.readline().split()[0])
        
        mesh_file.readline()

        global_node_id_type_zero = np.zeros((nen_type_zero,numel_type_zero),int)

        for kk in range(0,numel_type_zero):
            temp_string_list = mesh_file.readline().split()
            for pp in range(0, nen_type_zero):
                global_node_id_type_zero[pp,kk] = temp_string_list[pp]

        mesh_file.readline()
        mesh_file.readline()
        mesh_file.readline()

        for kk in range(0,numel_type_zero):
            vertex_entity_hash[kk+1] = int(mesh_file.readline().split()[0]) 

        elem_con_file_type_zero.write("type_of_element\n" + "0 -> point\n\n")
        elem_con_file_type_zero.write("type_of_entity\n" + "0 -> vertex\n\n")
        elem_con_file_type_zero.write("number_of_entities\n" + str(max(vertex_entity_hash.values())+1) + "\n\n")
        elem_con_file_type_zero.write("number_of_elements\n" + str(numel_type_zero) + "\n\n")
        elem_con_file_type_zero.write("number_of_nodes_per_element\n" + str(nen_type_zero) + "\n\n")
        elem_con_file_type_zero.write("entity_id  " +  "element_id  " + "global_node_ids  " + "x_node_i  " + "y_node_i  " + "z_node_i\n")

        for kk in range(0,numel_type_zero):       
            elem_con_file_type_zero.write(str(vertex_entity_hash[kk+1]) + "    ")
            elem_con_file_type_zero.write(str(kk+1) + "    ")

            for pp in range(0,nen_type_zero):
                elem_con_file_type_zero.write(str(global_node_id_type_zero[pp,kk]+1) + "  ")

            for pp in range(0,nen_type_zero):
                for gg in range(3):
                    elem_con_file_type_zero.write("  " + "%.5f" %(node_hash[global_node_id_type_zero[pp,kk]+1][gg]))
                elem_con_file_type_zero.write("        ")
            elem_con_file_type_zero.write('\n')

    if ("2 # number of nodes per element" in line):
        nen_type_one   = int(line.split()[0])
        numel_type_one = int(mesh_file.readline().split()[0])
        
        mesh_file.readline()

        global_node_id_type_one = np.zeros((nen_type_one,numel_type_one),int)

        for kk in range(0,numel_type_one):
            temp_string_list = mesh_file.readline().split()
            for pp in range(0, nen_type_one):
                global_node_id_type_one[pp,kk] = temp_string_list[pp]

        mesh_file.readline()
        mesh_file.readline()
        mesh_file.readline()

        for kk in range(0,numel_type_one):
            edge_entity_hash[kk+1] = int(mesh_file.readline().split()[0]) 

        elem_con_file_type_one.write("type_of_element\n" + "1 -> line\n\n")
        elem_con_file_type_one.write("type_of_entity\n" + "1 -> edge\n\n")
        elem_con_file_type_one.write("number_of_entities\n" + str(max(edge_entity_hash.values())+1) + "\n\n")
        elem_con_file_type_one.write("number_of_elements\n" + str(numel_type_one) + "\n\n")
        elem_con_file_type_one.write("number_of_nodes_per_element\n" + str(nen_type_one) + "\n\n")
        elem_con_file_type_one.write("entity_id  " + "element_id  " + "global_node_ids  " + "x_node_i  " + "y_node_i  " + "z_node_i\n")

        for kk in range(0,numel_type_one):
            elem_con_file_type_one.write(str(edge_entity_hash[kk+1]) + "    ")
            elem_con_file_type_one.write(str(kk+1) + "    ")

            for pp in range(0,nen_type_one):
                elem_con_file_type_one.write(str(global_node_id_type_one[pp,kk]+1) + "  ")

            for pp in range(0,nen_type_one):
                for gg in range(3):
                    elem_con_file_type_one.write("  " + "%.5f" %(node_hash[global_node_id_type_one[pp,kk]+1][gg]))
                elem_con_file_type_one.write("        ")
            elem_con_file_type_one.write('\n')
        
    if ("3 # number of nodes per element" in line):
        nen_type_two   = int(line.split()[0])
        numel_type_two = int(mesh_file.readline().split()[0])
        
        mesh_file.readline()

        global_node_id_type_two = np.zeros((nen_type_two,numel_type_two),int)

        for kk in range(0,numel_type_two):
            temp_string_list = mesh_file.readline().split()
            for pp in range(0, nen_type_two):
                global_node_id_type_two[pp,kk] = temp_string_list[pp]

        mesh_file.readline()
        mesh_file.readline()
        mesh_file.readline()

        for kk in range(0,numel_type_two):
            face_entity_hash[kk+1] = int(mesh_file.readline().split()[0]) 

        elem_con_file_type_two.write("type_of_element\n" + "2 -> triangle\n\n")
        elem_con_file_type_two.write("type_of_entity\n" + "2 -> face\n\n")
        elem_con_file_type_two.write("number_of_entities\n" + str(max(face_entity_hash.values())+1) + "\n\n")
        elem_con_file_type_two.write("number_of_elements\n" + str(numel_type_two) + "\n\n")
        elem_con_file_type_two.write("number_of_nodes_per_element\n" + str(nen_type_two) + "\n\n")
        elem_con_file_type_two.write("entity_id  " + "element_id  " + "global_node_ids  " + "x_node_i  " + "y_node_i  " + "z_node_i\n")

        for kk in range(0,numel_type_two):
            elem_con_file_type_two.write(str(face_entity_hash[kk+1]) + "    ")
            elem_con_file_type_two.write(str(kk+1) + "    ")

            for pp in range(0,nen_type_two):
                elem_con_file_type_two.write(str(global_node_id_type_two[pp,kk]+1) + "  ")

            for pp in range(0,nen_type_two):
                for gg in range(3):
                    elem_con_file_type_two.write("  " + "%.5f" %(node_hash[global_node_id_type_two[pp,kk]+1][gg]))
                elem_con_file_type_two.write("        ")
            elem_con_file_type_two.write('\n')

    if ("4 # number of nodes per element" in line):
        nen_type_three   = int(line.split()[0])
        numel_type_three = int(mesh_file.readline().split()[0])

        mesh_file.readline()

        global_node_id_type_three = np.zeros((nen_type_three,numel_type_three),int)

        for kk in range(0,numel_type_three):
            temp_string_list = mesh_file.readline().split()
            for pp in range(0, nen_type_three):
                global_node_id_type_three[pp,kk] = temp_string_list[pp]

        mesh_file.readline()
        mesh_file.readline()
        mesh_file.readline()

        for kk in range(0,numel_type_three):
            domain_entity_hash[kk+1] = int(mesh_file.readline().split()[0]) 

        elem_con_file_type_three.write("type_of_element\n" + "3 -> tetrahedron\n\n")
        elem_con_file_type_three.write("type_of_entity\n" + "3 -> domain\n\n")
        elem_con_file_type_three.write("number_of_entities\n" + str(max(domain_entity_hash.values())+1) + "\n\n")
        elem_con_file_type_three.write("number_of_elements\n" + str(numel_type_three) + "\n\n")
        elem_con_file_type_three.write("number_of_nodes_per_element\n" + str(nen_type_three) + "\n\n")
        elem_con_file_type_three.write("entity_id  " + "element_id  " + "global_node_ids  " + "x_node_i  " + "y_node_i  " + "z_node_i\n")

        for kk in range(0,numel_type_three):
            elem_con_file_type_three.write(str(domain_entity_hash[kk+1]) + "    ")
            elem_con_file_type_three.write(str(kk+1) + "    ")

            for pp in range(0,nen_type_three):
                elem_con_file_type_three.write(str(int(global_node_id_type_three[pp,kk])) + "  ")

            for pp in range(0,nen_type_three):
                for gg in range(3):
                    elem_con_file_type_three.write("  " + "%.5f" %(node_hash[global_node_id_type_three[pp,kk]+1][gg]))
                elem_con_file_type_three.write("        ")
            elem_con_file_type_three.write('\n')
        
mesh_file.close()
meshpoints_file.close()
elem_con_file_type_zero.close()
elem_con_file_type_one.close()
elem_con_file_type_two.close()
elem_con_file_type_three.close()

node_pair_xtox = {}
node_pair_ytoy = {}
node_pair_ztoz = {}

x_face1_elements = []
x_face2_elements = []

y_face1_elements = []
y_face2_elements = []

z_face1_elements = []
z_face2_elements = []

for key, value in face_entity_hash.items():
    if value == 0: x_face1_elements.append(key)
    if value == 5: x_face2_elements.append(key)
    if value == 1: y_face1_elements.append(key)
    if value == 4: y_face2_elements.append(key)
    if value == 3: z_face1_elements.append(key)
    if value == 2: z_face2_elements.append(key)

for ii in range(len(x_face1_elements)):
    xface1.write(str(x_face1_elements[ii]) + ' ')

for ii in range(len(x_face2_elements)):
    xface2.write(str(x_face2_elements[ii]) + ' ')

for ii in range(len(y_face1_elements)):
    yface1.write(str(y_face1_elements[ii]) + ' ')

for ii in range(len(y_face2_elements)):
    yface2.write(str(y_face2_elements[ii]) + ' ')

for ii in range(len(z_face1_elements)):
    zface1.write(str(z_face1_elements[ii]) + ' ')

for ii in range(len(z_face2_elements)):
    zface2.write(str(z_face2_elements[ii]) + ' ')

# Build xx pairing
for elem1 in x_face1_elements:
    for elem2 in x_face2_elements:
        for ii in range(3):
            node1 = global_node_id_type_two[ii,elem1-1] + 1
            for jj in range(3):
                node2 = global_node_id_type_two[jj,elem2-1] + 1
                if ((abs(node_hash[node1][1]-node_hash[node2][1])<tol)and(abs(node_hash[node1][2]-node_hash[node2][2])<tol)):
                    node_pair_xtox[node1] = node2
for key, value in node_pair_xtox.items():
    node_pairing_xx.write(str(key) + "  " + str(value) + '\n')

# Build yy pairing    
for elem1 in y_face1_elements:
    for elem2 in y_face2_elements:
        for ii in range(3):
            node1 = global_node_id_type_two[ii,elem1-1] + 1
            for jj in range(3):
                node2 = global_node_id_type_two[jj,elem2-1] + 1
                if ((abs(node_hash[node1][0]-node_hash[node2][0])<tol)and(abs(node_hash[node1][2]-node_hash[node2][2])<tol)):
                    node_pair_ytoy[node1] = node2
for key, value in node_pair_ytoy.items():
    node_pairing_yy.write(str(key) + "  " + str(value) + '\n')

# Build zz pairing
for elem1 in z_face1_elements:
    for elem2 in z_face2_elements:
        for ii in range(3):
            node1 = global_node_id_type_two[ii,elem1-1] + 1
            for jj in range(3):
                node2 = global_node_id_type_two[jj,elem2-1] + 1
                if ((abs(node_hash[node1][0]-node_hash[node2][0])<tol)and(abs(node_hash[node1][1]-node_hash[node2][1])<tol)):
                    node_pair_ztoz[node1] = node2
for key, value in node_pair_ztoz.items():
    node_pairing_zz.write(str(key) + "  " + str(value) + '\n')

exit()
