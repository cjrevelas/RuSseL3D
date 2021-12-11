import numpy as np
import math as m
import sys


TWO_PI  = np.pi * 2.0
PI      = np.pi
HALF_PI = np.pi * 0.5


class bias():
    def __init__(self, phi, the, mag, sig):
        self.phi = phi
        self.the = the
        self.mag = mag
        self.sig = sig


class nanoparticle():
    def __init__(self, r0, rad):
       self.r0   = r0
       self.rad  = rad
       self.gps  = []
       self.n_gp = 0


    def add_gp(self, phi, the):
        self.gps.append([phi, the])
        self.n_gp = len(self.gps)

        return


    def gp_pos(self, id):
        phi = self.gps[id][0] - HALF_PI
        the = self.gps[id][1] - PI

        rr    = [0.0, 0.0, 0.0]
        rr[0] = self.r0[0] + self.rad * np.sin(phi) * np.cos(the)
        rr[1] = self.r0[1] + self.rad * np.sin(phi) * np.sin(the)
        rr[2] = self.r0[2] + self.rad * np.cos(phi)

        return rr


def get_great_angle(phi1, the1, phi2, the2):
    dthe = abs(the2 - the1)
    ds   = np.arccos( np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(dthe) )

    return ds


def sample_random_point():
    the = np.random.uniform(0.0, TWO_PI) - PI
    v   = np.random.uniform(0.0, 1.0)
    phi = np.arccos(1.0 - 2.0*v) - HALF_PI

    return phi, the


def get_points_from_equidis(n_gp):
    alpha   = 4.0 * np.pi / n_gp
    da      = np.sqrt(alpha)
    m_theta = int(round(np.pi/d))
    d_theta = np.pi/m_theta
    d_phi   = alpha/d_theta
    coords  = []

    for ii in range(0, m_theta):
        theta = np.pi * (ii + 0.5)/m_theta
        m_phi = int(round(2 * np.pi * np.sin(theta)/d_phi))
        for jj in range(0, m_phi):
            phi = (2.0 * np.pi * float(jj))/float(m_phi)
            coords.append([theta-HALF_PI, phi-PI])

    return coords


#set parameter values (Angstrom units)
n_np   = 1
np_pos = [0.0, 0.0, 0.0]
radius = 40.0 + 4.4
n_gp   = 80

min_gp_dist = 12.0

n_the = 400 # grid points along the theta axis (-pi,+pi)
n_phi = 200 # grid points along the phi axis   (-pi/2, +pi/2)

# set the position, magnitude and deviation of the biasing poles: phi, the, mag, sig
biases = []

biases.append(bias(+HALF_PI, 0.0, -0.7, 1.0))
biases.append(bias(-HALF_PI, 0.0, -0.7, 1.0))

#coords = get_points_from_equidis(n_gp)
#for coord in coords:
#     print(coord)
#     [phi, the] = coord
#     biases.append(bias(phi, the, 0.99, 1.0))

prob_background = 1.0

# construct the number of phi/theta grid points
d_the    = TWO_PI / n_the
d_phi    = PI     / n_phi
grid_the = [d_the * ii - PI for ii in range(n_the)]
grid_phi = [d_phi * ii - HALF_PI for ii in range(n_phi)]

# calculate the probability map
prob_ref = [[prob_background for ii in range(n_the)] for jj in range(n_phi)]

for bias in biases:
    phi1 = bias.phi
    the1 = bias.the
    mag = bias.mag
    sig = bias.sig
    print(bias.phi, bias.the)

    for ii in range(n_phi):
        for jj in range(n_the):
            phi2 = grid_phi[ii]
            the2 = grid_the[jj]

            ds = get_great_angle(phi1, the1, phi2, the2)
            dr = radius * ds

            prob = mag * np.exp( -0.5 * pow(dr/sig,2) )

            prob_ref[ii][jj] += prob

#check whether the probability lies within the range [0,1] bounds
for ii in range(n_phi):
    for jj in range(n_the):
        if prob_ref[ii][jj] <= -0.0001:
            print("Prob for phi "+str(grid_phi[ii])+", the "+str(grid_the[jj])+" is outside the bound [0,1] ("+str(prob_ref[ii][jj])+")")
            prob_ref[ii][jj] = 0.0

        if prob_ref[ii][jj] >= 1.0001:
            print("Prob for phi "+str(grid_phi[ii])+", the "+str(grid_the[jj])+" is outside the bound [0,1] ("+str(prob_ref[ii][jj])+")")
            prob_ref[ii][jj] = 1.0

#export the probability map
prob_map_file = open("o.prob_map.dat",'w')

prob_map_file.write("# ")
for jj in range(n_the):
    prob_map_file.write(str(grid_the[jj]) + ' ')
prob_map_file.write('\n')

for ii in range(n_phi):
    prob_map_file.write(str(grid_phi[ii]) + ' ')
    for jj in range(n_the):
        #if norm != 0: prob_ref[ii][jj] = (prob_ref[ii][jj] - prob_min ) / norm
        prob_map_file.write(str(prob_ref[ii][jj]) + ' ')
    prob_map_file.write('\n')

prob_map_file.close()

#create an instance of the np class and start inserting the grafting points
np1 = nanoparticle(np_pos, radius)

n_gp_cur = 0


np1.add_gp(0,0)

while np1.n_gp < n_gp:
    [phi, the] = sample_random_point()

    iphi = int((phi+HALF_PI) / d_phi)
    ithe = int((the+PI)      / d_the)

    #print(the, grid_the[ithe], the-grid_the[ithe], phi, grid_phi[iphi], phi - grid_phi[iphi])

    prob_insert = prob_ref[iphi][ithe]

    if np.random.rand() < prob_insert:
        #check that gp distance is above minimum
        insert = True
        for gp in np1.gps:
            dist = get_great_angle(phi, the, gp[0], gp[1]) * np1.rad
            #print(dist, min_gp_dist)
            if dist < min_gp_dist:
                insert = False
                break
        if insert:
           np1.add_gp(phi, the)
           print(np1.n_gp)


#export the grafting points
gnodes_file = open("o.gnodes.dist.lammpstrj",'w')
gnodes_file.write("ITEM: TIMESTEP\n")
gnodes_file.write("0\n")
gnodes_file.write("ITEM: NUMBER OF ATOMS\n")
gnodes_file.write("%d\n" % (n_gp+1))
gnodes_file.write("ITEM: BOX BOUNDS pp pp pp\n")
gnodes_file.write("%f %f\n" % (np_pos[0]-1.2*radius, np_pos[0]+1.2*radius))
gnodes_file.write("%f %f\n" % (np_pos[1]-1.2*radius, np_pos[0]+1.2*radius))
gnodes_file.write("%f %f\n" % (np_pos[2]-1.2*radius, np_pos[0]+1.2*radius))
gnodes_file.write("ITEM: ATOMS id type x y z \n")

for id in range(np1.n_gp):
    rr = np1.gp_pos(id)
    gnodes_file.write("%d %d %f %f %f\n" % (id+1, 0, rr[0], rr[1], rr[2]))
gnodes_file.write("%d %d %f %f %f\n" % (id+2, 1, np_pos[0], np_pos[1], np_pos[2]))
gnodes_file.close()

gnodes_file = open("gpoints.txt",'w')
for id in range(np1.n_gp):
    rr = np1.gp_pos(id)
    gnodes_file.write("%f %f %f\n" % (rr[0], rr[1], rr[2]))
gnodes_file.close()
