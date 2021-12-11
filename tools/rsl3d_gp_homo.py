import numpy as np
import math as m
import sys

TWO_PI  = np.pi * 2.0
PI      = np.pi
HALF_PI = np.pi * 0.5

class nanoparticle():
   def __init__(self, r0, rad):
      self.r0  = r0
      self.rad = rad
      self.gps = []
      self.n_gp = 0
   def add_gp(self, phi, the):
      self.gps.append([phi, the])
      self.n_gp = len(self.gps)
   def gp_pos(self, id):
      phi = self.gps[id][0] - HALF_PI
      the = self.gps[id][1] - PI

      rr    = [0.0, 0.0, 0.0]
      rr[0] = self.r0[0] + self.rad * np.sin(phi) * np.cos(the)
      rr[1] = self.r0[1] + self.rad * np.sin(phi) * np.sin(the)
      rr[2] = self.r0[2] + self.rad * np.cos(phi)
      return rr

def get_points_from_equidis(n_gp):
   alpha = 4 * np.pi / n_gp
   d = np.sqrt(alpha)
   m_theta = int(round(np.pi/d))
   d_theta = np.pi/m_theta
   d_phi   = alpha/d_theta
   coords = []
   for ii in range(0, m_theta):
       theta = np.pi * (ii + 0.5)/m_theta
       m_phi = int(round(2 * np.pi * np.sin(theta)/d_phi))
       for jj in range(0, m_phi):
           phi = (2 * np.pi * jj)/m_phi
           coords.append([theta-HALF_PI, phi-PI])
   return coords

# set the radious of the NP
n_np     = 1    # number of nanoparticles
np_pos   = [0.0, 0.0, 0.0]
radious  = 80+4.4   # radii of the nanoparticle
n_gp     = 500  # number of grafting points

n_the    = 600  # grid points along the theta axis (-pi,+pi)
n_phi    = 300 # grid points along the phi axis   (-pi/2, +pi/2)

prob_background = 0.0

#construct the number of phi/theta grid points
d_the    = TWO_PI / n_the
d_phi    = PI     / n_phi
grid_the = [d_the * ii - PI for ii in range(n_the)]
grid_phi = [d_phi * ii - HALF_PI for ii in range(n_phi)]

#calculate the probability map
prob_ref = [[prob_background for ii in range(n_the)] for jj in range(n_phi)]

# initialize the np class
np1 = nanoparticle(np_pos, radious)
# start the insertion of the grafting points

coords = get_points_from_equidis(n_gp)
for coord in coords:
   [phi, the] = coord
   print(phi, the)
   iphi = int( (phi+HALF_PI) / d_phi)
   ithe = int( (the+PI)      / d_the)
   prob_ref[iphi][ithe] += 1
   np1.add_gp(phi, the)

# output the grafting points
foo = open("o.gnodes.dist.lammpstrj","w")
foo.write("ITEM: TIMESTEP\n")
foo.write("0\n")
foo.write("ITEM: NUMBER OF ATOMS\n")
foo.write("%d\n" % (n_gp+1))
foo.write("ITEM: BOX BOUNDS pp pp pp\n")
foo.write("%f %f\n" % (np_pos[0]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("%f %f\n" % (np_pos[1]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("%f %f\n" % (np_pos[2]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("ITEM: ATOMS id type x y z \n")
for id in range(np1.n_gp):
   rr = np1.gp_pos(id)
   foo.write("%d %d %f %f %f\n" % (id+1, 0, rr[0], rr[1], rr[2]))
   id += 1
foo.write("%d %d %f %f %f\n" % (id+1, 1, np_pos[0], np_pos[1], np_pos[2]))
foo.close()

# output the probability map
foo = open("o.prob_map.dat","w")

foo.write("# ")
for jj in range(n_the):
    foo.write(str(grid_the[jj])+" ")
foo.write("\n")

for ii in range(n_phi):
   foo.write(str(grid_phi[ii])+" ")
   for jj in range(n_the):
#      if norm != 0: prob_ref[ii][jj] = (prob_ref[ii][jj] - prob_min ) / norm
      foo.write(str(prob_ref[ii][jj])+" ")
   foo.write("\n")

foo.close()
