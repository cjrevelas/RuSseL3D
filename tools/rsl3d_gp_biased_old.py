import numpy as np
import math as m
import sys

TWO_PI = 2.0*np.pi
PI = np.pi
HALF_PI = np.pi*0.5

class bias():
   def __init__(self, phi, the, mag, sig):
      self.phi = phi
      self.the = the
      self.mag = mag
      self.sig = sig

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
      phi = self.gps[id][0]
      the = self.gps[id][1]
      rr    = [0.0, 0.0, 0.0]
      rr[0] = self.r0[0] + self.rad * np.sin(phi) * np.cos(the)
      rr[1] = self.r0[1] + self.rad * np.sin(phi) * np.sin(the)
      rr[2] = self.r0[2] + self.rad * np.cos(phi)
      return rr

def sample_random_point():
   the = np.random.uniform(0,TWO_PI);
   v = np.random.uniform(0,1);
   phi = np.arccos(1-2*v);
   return phi, the

####################################################
# Parameter section
#
# set the radious of the NP
n_np     = 1    # number of nanoparticles
np_pos   = [10.0, 10.0, 10.0]
radious  = 10   # radii of the nanoparticle
n_gp     = 500  # number of grafting points

n_the    = 100 # grid points along the theta axis (-pi,+pi)
n_phi    = 50  # grid points along the phi axis   (-pi/2, +pi/2)

# set the position, magnitude and deviation of the biasing poles
#                  phi      ,the      ,mag      , sig
biases = []
#biases.append(bias(HALF_PI  ,HALF_PI  ,1        ,   0.5))
biases.append(bias(HALF_PI  ,0.0      , 1.0         ,   5))
biases.append(bias(-HALF_PI  ,0.0      ,0.0         ,   5))

prob_background = 0.05
####################################################
#
# construct the number of phi/theta grid points
d_the    = TWO_PI / n_the
d_phi    = PI     / n_phi
grid_the = [d_the * ii for ii in range(n_the)]
grid_phi = [d_phi * ii for ii in range(n_phi)]

# calculate the probability map
prob_ref = [[prob_background for ii in range(n_the)] for jj in range(n_phi)]

for bias in biases:
   phi1 = bias.phi
   the1 = bias.the
   mag = bias.mag-prob_background
   sig = bias.sig-prob_background

   for ii in range(n_phi):
      for jj in range(n_the):
         phi2 = grid_phi[ii]
         the2 = grid_the[jj]

         def get_great_angle(phi1, the1, phi2, the2):
            dphi = abs(phi2 - phi1)
            ds = np.arccos( np.sin(the1) * np.sin(the2) + np.cos(the1) * np.cos(the2) * np.cos(dphi) )
            return ds
         ds = get_great_angle(phi1, the1, phi2, the2)
         #dthe = abs(the2 - the1)
         #ds = np.arccos( np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(dthe) )
         dr = radious * ds
         prob = mag * np.exp( -0.5 * pow((dr)/sig,2) )

         prob_ref[ii][jj] += prob

# check whether the probability lies within the [0,1] bounds
for ii in range(n_phi):
   for jj in range(n_the):
      if prob_ref[ii][jj] < 0 or prob_ref[ii][jj] > 1:
         print("Prob for phi "+str(grid_phi[ii])+", the "+str(grid_the[jj])+" is outside the bound [0,1] ("+str(prob_ref[ii][jj])+")")
         sys.exit()
#flattened = [val for sublist in prob_ref for val in sublist]
#
#prob_max = max(flattened)
#prob_min = min(flattened)
#
#norm = prob_max - prob_min

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

# initialize the np class

np1 = nanoparticle(np_pos, radious)
# start the insertion of the grafting points

n_gp_cur = 0
while n_gp_cur < n_gp:
   [phi, the] = sample_random_point()
   
   iphi = int(phi/d_phi)
   ithe = int(the/d_the)

   prob_insert = prob_ref[iphi][ithe]

   if np.random.rand() < prob_insert:
      np1.add_gp(phi, the)
      n_gp_cur += 1

# output the grafting points
foo = open("o.gnodes.dist.lammpstrj","w")
foo.write("ITEM: TIMESTEP\n")
foo.write("0\n")
foo.write("ITEM: NUMBER OF ATOMS\n")
foo.write("%d\n" % (n_gp))
foo.write("ITEM: BOX BOUNDS pp pp pp\n")
foo.write("%f %f\n" % (np_pos[0]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("%f %f\n" % (np_pos[1]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("%f %f\n" % (np_pos[2]-1.2*radious, np_pos[0]+1.2*radious))
foo.write("ITEM: ATOMS id type x y z \n")
id = 1
for id in range(np1.n_gp):
   rr = np1.gp_pos(id)
   foo.write("%d %d %f %f %f\n" % (id, 0, rr[0], rr[1], rr[2]))
   id += 1
foo.close()
