import numpy as np
import os

coord_file = "coord_nodes"
el_file = "el_nodes"
central_node = 6146

def tri_area(a, b, c) :
    return 0.5 * np.linalg.norm( np.cross( b-a, c-a ) )

def tetvol(a, b, c, d):
    return (abs(np.linalg.det((a-b, b-c, c-d))) / 6.0)

#test tetvol
#a = np.array([0.0, 0.0, 0.0])
#b = np.array([0.0, 0.0, 2.0])
#c = np.array([0.0, 2.0, 0.0])
#d = np.array([2.0, 0.0, 0.0])
#print(tetvol(a, b, c, d)) # = 1.3333..."

#test tri_area
#a = np.array([0.0, 0.0, 0.0])
#b = np.array([2.0, 0.0, 0.0])
#c = np.array([0.0, 2.0, 0.0])
#print(tri_area(a, b, c)) # = 2.000...#
#quit()

#
# Read the file with the coordinates
#
foo = open(coord_file,"r")
# read number of nodes
line = foo.readline().split()
n_nodes = int(line[0])
print("number of nodes:", n_nodes)
# skip header
line = foo.readline()
# read coordinates
node_coord = {}
for ii in range(n_nodes):
   line = foo.readline().split()
   id = int(line[0])
   xx = float(line[1])
   yy = float(line[2])
   zz = float(line[3])
   node_coord[id] = np.array([xx, yy, zz])
foo.close()
#
# Read the file with the nodes per element
#
foo = open(el_file,"r")
# read number of elements
line = foo.readline().split()
n_el = int(line[0])
print("number of elements:", n_el)
# skip header
line = foo.readline()
# read coordinates
el_nodes = {}
for ii in range(n_el):
   line = foo.readline().split()
   id = int(line[0])
   ii = int(line[1])
   jj = int(line[2])
   kk = int(line[3])
   ll = int(line[4])
   el_nodes[id] = [ii, jj, kk, ll]
foo.close()
#
# Compute the volume of each element, the areas of faces
# and the contributing volume to each node
#
vol_central_node = 0.0
print("element ii jj kk ll Ei Ej Ek El Vi Vj Vk Vl vol vol_node")
el_vol = {}
for el in el_nodes:
   ii = el_nodes[el][0]
   jj = el_nodes[el][1]
   kk = el_nodes[el][2]
   ll = el_nodes[el][3]

   rii = node_coord[ii]
   rjj = node_coord[jj]
   rkk = node_coord[kk]
   rll = node_coord[ll]

   vol = tetvol(rii, rjj, rkk, rll)

   area = {}
   area[ii] = tri_area(rjj, rkk, rll)
   area[jj] = tri_area(rkk, rll, rii)
   area[kk] = tri_area(rll, rii, rjj)
   area[ll] = tri_area(rii, rjj, rkk)

   tot_area = area[ii] + area[jj] + area[kk] + area[ll]

   vol_node = {}
   vol_node[ii] = vol * area[ii] / tot_area
   vol_node[jj] = vol * area[jj] / tot_area
   vol_node[kk] = vol * area[kk] / tot_area
   vol_node[ll] = vol * area[ll] / tot_area

   vol_central_node += vol_node[central_node]

   print(el, ii, jj, kk, ll, \
         area[ii], area[jj], area[kk], area[ll], \
         vol_node[ii], vol_node[jj], vol_node[kk], vol_node[ll], \
         vol, vol_node[central_node])

print("\nnode volume: ", vol_central_node)
