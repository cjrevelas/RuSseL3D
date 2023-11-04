import numpy as np
import sys

kHighVal     = 1e10
kTol         = 1e-8
adjust_nodes = True
periodicity  = [True, True, True]
path_mesh    = "in.mesh.original"

# name of the adjusted mesh
path_mesh_new = "in.mesh.refined"

# open the input file
foo = open(path_mesh, 'r')

n_node = 0
lines = []
while True:
   line = foo.readline()
   lines.append(line)
   if '# number of mesh points' in line:
      n_node = int(line.split()[0])
   if '# Mesh point coordinates' in line:
      break

nodes = []
box_lo = [+kHighVal, +kHighVal, +kHighVal]
box_hi = [-kHighVal, -kHighVal, -kHighVal]
for inode in range(n_node):
   line_split = foo.readline().split()
   r0 = float(line_split[0])
   r1 = float(line_split[1])
   r2 = float(line_split[2])

   # find the box bounds
   box_lo[0] = min(r0, box_lo[0])
   box_lo[1] = min(r1, box_lo[1])
   box_lo[2] = min(r2, box_lo[2])
   box_hi[0] = max(r0, box_hi[0])
   box_hi[1] = max(r1, box_hi[1])
   box_hi[2] = max(r2, box_hi[2])

   # store the node
   nodes.append([r0, r1, r2])

nodes_edge = [ [[],[]], [[],[]], [[],[]]] # [dim][lo/hi][id]

id = 0
for [r0, r1, r2] in nodes:
   if abs(r0 - box_lo[0]) < kTol: nodes_edge[0][0].append(id)
   if abs(r1 - box_lo[1]) < kTol: nodes_edge[1][0].append(id)
   if abs(r2 - box_lo[2]) < kTol: nodes_edge[2][0].append(id)
   if abs(r0 - box_hi[0]) < kTol: nodes_edge[0][1].append(id)
   if abs(r1 - box_hi[1]) < kTol: nodes_edge[1][1].append(id)
   if abs(r2 - box_hi[2]) < kTol: nodes_edge[2][1].append(id)
   id += 1

print("n_node  ", n_node)
print(box_lo[0], box_hi[0], len(nodes_edge[0][0]), len(nodes_edge[0][1]), ' # xlo xhi nlo nhi')
print(box_lo[1], box_hi[1], len(nodes_edge[1][0]), len(nodes_edge[1][1]), ' # ylo yhi nlo nhi')
print(box_lo[2], box_hi[2], len(nodes_edge[2][0]), len(nodes_edge[2][1]), ' # zlo zhi nlo nhi')

# Check whether the mumber of nodes at each face is the same
if len(nodes_edge[0][0]) != len(nodes_edge[0][1]):
  print("Error: difference in number of nodes along axis 0")
if len(nodes_edge[1][0]) != len(nodes_edge[1][1]):
  print("Error: difference in number of nodes along axis 1")
if len(nodes_edge[2][0]) != len(nodes_edge[2][1]):
  print("Error: difference in number of nodes along axis 2")


def CheckPeriodicity(nodes, nodes_edge, m0, adjust_nodes):

   # dimensions orthogonal to m0
   m1 = (m0+1)%3
   m2 = (m0+2)%3

   max_dij = -kHighVal
   for ii_base in nodes_edge[m0][0]:
      dij     = kHighVal
      ii_pair = 0
      r1      = nodes[ii_base][m1]
      r2      = nodes[ii_base][m2]
      for jj in nodes_edge[m0][1]:
         d1  = r1 - nodes[jj][m1]
         d2  = r2 - nodes[jj][m2]
         d12 = np.sqrt(d1*d1 + d2*d2)
         if d12 < dij:
            ii_pair = jj
            dij     = d12
      max_dij = max(max_dij, dij)
      if dij > kTol:
         print("   pair mismatch: ",ii_base, ii_pair, dij)
         #print("   pair mismatch: ",ii_base, ii_pair, dij, nodes[ii], nodes[ii_pair])
         if adjust_nodes:
            r1_av = 0.5*(nodes[ii_base][m1] + nodes[ii_pair][m1])
            r2_av = 0.5*(nodes[ii_base][m2] + nodes[ii_pair][m2])
            nodes[ii_base][m1] = r1_av
            nodes[ii_pair][m1] = r1_av
            nodes[ii_base][m2] = r2_av
            nodes[ii_pair][m2] = r2_av
   return max_dij

# Execute the periodicity check and node readjustement
for dim in range(3):
   if periodicity[dim]:
      print("axis:", dim)
      error = CheckPeriodicity(nodes, nodes_edge, dim, adjust_nodes)
      print("max_error: ", error)

# append adjusted coordinates
for node in nodes:
   line = "%.15f %.15f %.15f\n" %(node[0], node[1], node[2])
   lines.append(line)

# read the remaining lines of the original mesh file
while True:
   line = foo.readline()
   if not line:
      break
   lines.append(line)
foo.close()

# generate a the new mesh
goo = open(path_mesh_new, 'w')
for line in lines:
   goo.write(line)
goo.close()
