import numpy as np
import math as m
import sys

np.random.seed(seed=1234)

TWO_PI  = np.pi * 2.0
PI      = np.pi
HALF_PI = np.pi * 0.5

class Bias():
    def __init__(self, phi, the, mag, sig):
        self.phi = phi
        self.the = the
        self.mag = mag
        self.sig = sig

class Gp():
    def __init__(self, phi, the, rr):
       self.phi = phi
       self.the = the
       self.rr  = rr

class Nanoparticle():
    def __init__(self, r0, rad, rad_actual, n_gp_target):
       self.r0          = r0
       self.rad         = rad
       self.rad_actual  = rad_actual
       self.n_gp_target = n_gp_target
       self.gps         = []

    def AddGp(self, phi, the):
        phi_aux = phi - HALF_PI
        the_aux = the - PI
        rr  = []
        rr.append(self.r0[0] + self.rad * np.sin(phi_aux) * np.cos(the_aux))
        rr.append(self.r0[1] + self.rad * np.sin(phi_aux) * np.sin(the_aux))
        rr.append(self.r0[2] + self.rad * np.cos(phi_aux))
        self.gps.append(Gp(phi, the, rr))
        return

def GetGreatAngle(phi1, the1, phi2, the2): #arclength
    dthe = abs(the2 - the1)
    ds   = np.arccos( np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(dthe) )
    return ds

def SampleRandomPoint():
    the = np.random.uniform(0,TWO_PI) - PI
    v   = np.random.uniform(0,1)
    phi = np.arccos(1-2*v) - HALF_PI
    return phi, the

def GetPointsFromEquidist(n_gp):
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

def AssemblyProbabilityMap(radius, n_the, n_phi, grid_the, grid_phi, biases, prob_bgd):
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
                ds = GetGreatAngle(phi1, the1, phi2, the2)
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

    # Export the probability map
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

def WrapInBox(rr, bounds):
    box = [bounds[i][1] - bounds[i][0] for i in range(3)]
    rr_wrap = [rr[0] - bounds[0][0], rr[1] - bounds[1][0], rr[2] - bounds[2][0]]
    rr_wrap[0] -= box[0]*m.floor(rr_wrap[0]/box[0])
    rr_wrap[1] -= box[1]*m.floor(rr_wrap[1]/box[1])
    rr_wrap[2] -= box[2]*m.floor(rr_wrap[2]/box[2])
    rr_wrap[0] += bounds[0][0]
    rr_wrap[1] += bounds[1][0]
    rr_wrap[2] += bounds[2][0]
    return rr_wrap

#
# Interface
#

# Box bounds
Lbox   = 200.0
bounds = [[-Lbox/2, Lbox/2], [-Lbox/2, Lbox/2],[-Lbox/2, Lbox/2]] # [[xlo,xhi],[ylo,yhi],[zlo,zhi]]
box    = [bounds[i][1] - bounds[i][0] for i in range(3)]

wrap_gps     = True
vmd_friendly = False

# Grid properties
n_the    = 400
n_phi    = 200
d_the    = TWO_PI / n_the #(-pi,+pi)
d_phi    = PI     / n_phi #(-pi/2,+pi/2)
grid_the = [d_the * ii - PI for ii in range(n_the)]
grid_phi = [d_phi * ii - HALF_PI for ii in range(n_phi)]

# Same radious and number of grafted points
r_np_actual = 40.0
h_wall      = 4.0
h_gp        = 0.4
radius      = r_np_actual+h_wall+h_gp
radius_act  = r_np_actual
n_gp        = 30
min_gp_dist = 9.8

#
# BCC lattice
#
# Generate the nanoparticles
nanoparticles = []
nanoparticles.append(Nanoparticle([ 0.0         , 0.0         , 0.0         ], radius, radius_act, n_gp))
nanoparticles.append(Nanoparticle([-bounds[0][0],-bounds[1][0],-bounds[2][0]], radius, radius_act, n_gp))

# Generate the images of the nanoparticles
nanoparticle_images = []
nanoparticle_images.append(Nanoparticle([-bounds[0][0],-bounds[1][0],+bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([-bounds[0][0],+bounds[1][0],-bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([-bounds[0][0],+bounds[1][0],+bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([+bounds[0][0],-bounds[1][0],-bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([+bounds[0][0],-bounds[1][0],+bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([+bounds[0][0],+bounds[1][0],-bounds[2][0]], radius, radius_act, 0))
nanoparticle_images.append(Nanoparticle([+bounds[0][0],+bounds[1][0],+bounds[2][0]], radius, radius_act, 0))

# Generate the probability maps
# Same probability map for each nanoparticle
prob_bgd  = 1.0
prob_maps = []
for nanop in nanoparticles:
    prob_maps.append(AssemblyProbabilityMap(nanop.rad, n_the, n_phi, grid_the, grid_phi, [], prob_bgd))

# Grafting section
for ii in range(len(nanoparticles)):
    nanop    = nanoparticles[ii]
    prob_map = prob_maps[ii]
    while len(nanop.gps) < nanop.n_gp_target:
        [phi, the] = SampleRandomPoint()
        iphi = int((phi+HALF_PI) / d_phi)
        ithe = int((the+PI)      / d_the)
        prob_insert = prob_map[iphi][ithe]
        if np.random.rand() < prob_insert:
            insert = True
            for gp in nanop.gps:
                dist = GetGreatAngle(phi, the, gp.phi, gp.the) * nanop.rad_actual
                #print(dist, min_gp_dist)
                if dist < min_gp_dist:
                    insert = False
                    break
            if insert:
                nanop.AddGp(phi, the)
                print(len(nanop.gps), phi, the)

# total number of grafting points
n_gp_tot = 0
for nanop in nanoparticles:
    n_gp_tot += len(nanop.gps)

# Export a lammpstrj file
shift = [0.0, 0.0, 0.0]
if vmd_friendly:
  shift = [-bounds[0][0], -bounds[1][0], -bounds[2][0]]

gnodes_file    = open("o.gnodes.lammpstrj",'w')
gnodes_file.write("ITEM: TIMESTEP\n")
gnodes_file.write("0\n")
gnodes_file.write("ITEM: NUMBER OF ATOMS\n")
gnodes_file.write("%d\n" % (n_gp_tot+len(nanoparticles)+len(nanoparticle_images)))
gnodes_file.write("ITEM: BOX BOUNDS pp pp pp\n")
gnodes_file.write("%f %f\n" % (bounds[0][0] + shift[0], bounds[0][1] + shift[0]))
gnodes_file.write("%f %f\n" % (bounds[1][0] + shift[1], bounds[1][1] + shift[1]))
gnodes_file.write("%f %f\n" % (bounds[2][0] + shift[2], bounds[2][1] + shift[2]))
gnodes_file.write("ITEM: ATOMS id type x y z \n")
global_id = 0
for nanop in nanoparticles:
    for gp in nanop.gps:
        global_id += 1
        rr = gp.rr
        if wrap_gps: rr = WrapInBox(rr, bounds)
        gnodes_file.write("%d %d %f %f %f\n" % (global_id, 0, rr[0] + shift[0], rr[1] + shift[1], rr[2] + shift[2]))
for nanop in nanoparticles:
    global_id += 1
    rr = nanop.r0
    gnodes_file.write("%d %d %f %f %f\n" % (global_id, 1, rr[0] + shift[0], rr[1] + shift[1], rr[2] + shift[2]))
for nanop in nanoparticle_images:
    global_id += 1
    rr = nanop.r0
    gnodes_file.write("%d %d %f %f %f\n" % (global_id, 1, rr[0] + shift[0], rr[1] + shift[1], rr[2] + shift[2]))
gnodes_file.close()

# Export a gpoint file
gp_coords_file = open("o.gpoints",'w')
for nanop in nanoparticles:
    for gp in nanop.gps:
        rr = gp.rr
        if wrap_gps: rr = WrapInBox(rr, bounds)
        gp_coords_file.write("%f %f %f\n" % (rr[0], rr[1], rr[2]))
gp_coords_file.close()

# Export a mesh input file
mesh_file = open("o.mesh_input",'w')
mesh_file.write("BOX\n")
mesh_file.write("%f %f xlo xhi\n" % (bounds[0][0], bounds[0][1]))
mesh_file.write("%f %f ylo yhi\n" % (bounds[1][0], bounds[1][1]))
mesh_file.write("%f %f zlo zhi\n" % (bounds[2][0], bounds[2][1]))
mesh_file.write("\n")
mesh_file.write("N_NP\n")
mesh_file.write("%d\n" %(len(nanoparticles)+len(nanoparticle_images)))
mesh_file.write("X_CENTER Y_CENTER Z_CENTER\n")
for nanop in nanoparticles:
    mesh_file.write("%f %f %f\n" % (nanop.r0[0], nanop.r0[1], nanop.r0[2]))
for nanop in nanoparticle_images:
    mesh_file.write("%f %f %f\n" % (nanop.r0[0], nanop.r0[1], nanop.r0[2]))
mesh_file.write("\n")
mesh_file.write("N_GP\n")
mesh_file.write("%d\n" %(n_gp_tot))
mesh_file.write("X_GP Y_GP Z_GP\n")
for nanop in nanoparticles:
    for gp in nanop.gps:
        rr = gp.rr
        if wrap_gps: rr = WrapInBox(rr, bounds)
        mesh_file.write("%f %f %f\n" % (rr[0], rr[1], rr[2]))
mesh_file.close()


#
# Miscellaneus settings
#
# first nanoparticle  -> vertical gp configuration
#biases        = []
#biases.append(Bias(-HALF_PI, 0.0, +1.0, 5.0))
#biases.append(Bias(+HALF_PI, 0.0, +1.0, 5.0))
#prob_maps.append(AssemblyProbabilityMap(radius[0], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))
# second nanoparticle -> horizontal gp configuration
#biases        = []
#biases.append(Bias(0.0, 0.0, +1.0, 5.0))
#biases.append(Bias(PI, 0.0, +1.0, 5.0))
#prob_maps.append(AssemblyProbabilityMap(radius[1], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))
# Vertical (phi, theta)
#prob_bgd = 0.0
#biases = []
#biases.append(Bias(-HALF_PI, 0.0, +1.0, 5.0))
#biases.append(Bias(+HALF_PI, 0.0, +1.0, 5.0))
#prob_maps.append(AssemblyProbabilityMap(radius[0], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))
#prob_maps.append(AssemblyProbabilityMap(radius[1], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))
# Horizontal (phi, theta)
#prob_bgd = 0.0
#biases = []
#biases.append(Bias(0.0, 0.0, +1.0, 5.0))
#biases.append(Bias(PI, 0.0, +1.0, 5.0)) -> equivalent to (0, +-pi, +1.0, 5.0)
#prob_maps.append(AssemblyProbabilityMap(radius[0], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))
#prob_maps.append(AssemblyProbabilityMap(radius[1], n_the, n_phi, grid_the, grid_phi, biases, prob_bgd))

