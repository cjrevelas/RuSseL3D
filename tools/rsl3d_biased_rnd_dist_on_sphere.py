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

class nanop():
    def __init__(self, r0, rad, rad_actual):
       self.r0         = r0
       self.rad        = rad
       self.rad_actual = rad_actual
       self.gps        = []
       self.n_gp       = 0

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


def get_great_angle(phi1, the1, phi2, the2): #arclength
    dthe = abs(the2 - the1)
    ds   = np.arccos( np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(dthe) )

    return ds


def sample_random_point():
    the = np.random.uniform(0,TWO_PI) - PI;
    v   = np.random.uniform(0,1);
    phi = np.arccos(1-2*v) - HALF_PI;

    return phi, the


def get_points_from_equidis(n_gp):
    alpha   = 4 * np.pi / n_gp
    da      = np.sqrt(alpha)
    m_theta = int(round(np.pi/d))
    d_theta = np.pi/m_theta
    d_phi   = alpha/d_theta
    coords  = []

    for ii in range(0, m_theta):
        theta = np.pi * (ii + 0.5)/m_theta
        m_phi = int(round(2 * np.pi * np.sin(theta)/d_phi))
        for jj in range(0, m_phi):
            phi = (2 * np.pi * jj)/m_phi
            coords.append([theta-HALF_PI, phi-PI])

    return coords


def assembly_probability_map(radius, n_the, n_phi, grid_the, grid_phi, biases, prob_bgd):
    prob_ref = [[prob_bgd for ii in range(n_the)] for jj in range(n_phi)]

    for bias in biases:
        phi1 = bias.phi
        the1 = bias.the
        mag  = bias.mag
        sig  = bias.sig

        for ii in range(n_phi):
            for jj in range(n_the):
                phi2 = grid_phi[ii]
                the2 = grid_the[jj]

                ds = get_great_angle(phi1, the1, phi2, the2)
                dr = radius * ds

                prob = mag * np.exp(-0.5 * pow(dr/sig,2))

                prob_ref[ii][jj] += prob

    for ii in range(n_phi):
        for jj in range(n_the):
            if prob_ref[ii][jj] <= -0.0001:
                print("Prob for phi "+str(grid_phi[ii])+", the "+str(grid_the[jj])+" is outside the bound [0,1] ("+str(prob_ref[ii][jj])+")")
                prob_ref[ii][jj] = 0.0

            if prob_ref[ii][jj] >= 1.0001:
                print("Prob for phi "+str(grid_phi[ii])+", the "+str(grid_the[jj])+" is outside the bound [0,1] ("+str(prob_ref[ii][jj])+")")
                prob_ref[ii][jj] = 1.0

    # Normalize the probability map
    max_prob = 0.0
    for ii in range(n_phi):
        for jj in range(n_the):
            max_prob = max(max_prob, prob_ref[ii][jj])
    for ii in range(n_phi):
        for jj in range(n_the):
            prob_ref[ii][jj] /= max_prob

    prob_map_file = open("o.prob_map_" + str(prob_bgd),'w')

    prob_map_file.write("# ")
    for jj in range(n_the):
        prob_map_file.write(str(grid_the[jj]) + ' ')
    prob_map_file.write('\n')

    for ii in range(n_phi):
        prob_map_file.write(str(grid_phi[ii]) + ' ')
        for jj in range(n_the):
            prob_map_file.write(str(prob_ref[ii][jj]) + ' ')
        prob_map_file.write('\n')

    prob_map_file.close()

    return prob_ref

r_np_actual_1 = 20.0
r_np_actual_2 = 20.0
n_np          = 2  
h_wall        = 4.0
h_gp          = 0.4
radius        = [r_np_actual_1+h_wall+h_gp, r_np_actual_2+h_wall+h_gp]
radius_act    = [r_np_actual_1, r_np_actual_2]
positions     = [[-31.0,0.0,0.0],[+31.0,0.0,0.0]]
n_gp          = [15, 15]
min_gp_dist   = 9.8
n_the         = 400 
n_phi         = 200 
nanop_array   = []
biases        = []
maps          = []

d_the         = TWO_PI / n_the #(-pi,+pi)
d_phi         = PI     / n_phi #(-pi/2,+pi/2)
grid_the      = [d_the * ii - PI for ii in range(n_the)]
grid_phi      = [d_phi * ii - HALF_PI for ii in range(n_phi)]

prob_bgd = 0.0
biases.append(bias(-0.0, 0.0, +1.0, 5.0))
biases.append(bias(+PI, 0.0, +1.0, 5.0))
prob_ref1 = assembly_probability_map(radius[0], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd)
prob_ref2 = assembly_probability_map(radius[1], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd)

#prob_bgd = 1.0
#biases.append(bias(-HALF_PI, 0.0, -1.0, 32.0))
#biases.append(bias(+HALF_PI, 0.0, -1.0, 32.0))
#prob_ref1 = assembly_probability_map(radius[0], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd)

#prob_bgd = 0.9
#prob_ref  = assembly_probability_map(radius[1], n_the, n_phi, grid_the, grid_phi, [], prob_bgd)

maps.append(prob_ref1)
maps.append(prob_ref2)

for ii in range(len(positions)):
    nanop_array.append(nanop(positions[ii],radius[ii],radius_act[ii]))

for ii in range(len(nanop_array)):
    nanop_array[ii].add_gp(0,0)

for ii in range(len(nanop_array)):
    while nanop_array[ii].n_gp < n_gp[ii]:
        [phi, the] = sample_random_point()

        iphi = int((phi+HALF_PI) / d_phi)
        ithe = int((the+PI)      / d_the)

        prob_insert = maps[ii][iphi][ithe]

        if np.random.rand() < prob_insert:
            insert = True
            for gp in nanop_array[ii].gps:
                dist = get_great_angle(phi, the, gp[0], gp[1]) * nanop_array[ii].rad_actual #check the arguments of the get_great_angle function
                #print(dist, min_gp_dist)
                if dist < min_gp_dist:
                    insert = False
                    break
            if insert:
                nanop_array[ii].add_gp(phi, the)
                print(nanop_array[ii].n_gp)

n_gp_tot = 0
for ii in range(len(nanop_array)):
    n_gp_tot += n_gp[ii]

# Export the grafting points
gnodes_file    = open("o.gnodes.lammpstrj",'w')
gp_coords_file = open("o.gpoints",'w')
gnodes_file.write("ITEM: TIMESTEP\n")
gnodes_file.write("0\n")
gnodes_file.write("ITEM: NUMBER OF ATOMS\n")
gnodes_file.write("%d\n" % (n_gp_tot+len(nanop_array)))
gnodes_file.write("ITEM: BOX BOUNDS pp pp pp\n")
gnodes_file.write("%f %f\n" % (nanop_array[0].r0[0]-1.2*radius[0], nanop_array[1].r0[0]+1.2*radius[1]))
gnodes_file.write("%f %f\n" % (nanop_array[0].r0[1]-1.2*radius[0], nanop_array[1].r0[1]+1.2*radius[1]))
gnodes_file.write("%f %f\n" % (nanop_array[0].r0[2]-1.2*radius[0], nanop_array[1].r0[2]+1.2*radius[1]))
gnodes_file.write("ITEM: ATOMS id type x y z \n")

global_gp_id = 0
for nanop_id in range(len(nanop_array)):
    for id in range(nanop_array[nanop_id].n_gp):
        global_gp_id += 1 
        rr = nanop_array[nanop_id].gp_pos(id)
        gnodes_file.write("%d %d %f %f %f\n" % (global_gp_id, 0, rr[0], rr[1], rr[2]))
        gp_coords_file.write("%f %f %f\n" % (rr[0], rr[1], rr[2]))

for nanop_id in range(len(nanop_array)):
    gnodes_file.write("%d %d %f %f %f\n" % (nanop_id+global_gp_id+1, 1, nanop_array[nanop_id].r0[0], nanop_array[nanop_id].r0[1], nanop_array[nanop_id].r0[2]))
gnodes_file.close()
gp_coords_file.close()
