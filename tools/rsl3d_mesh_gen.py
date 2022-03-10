#----------------------------------------------------------------------------------------------------------------------------------#
import numpy as np
import os
import sys
#----------------------------------------------------------------------------------------------------------------------------------#
class Mesh:
    def __init__(self, planar, grafting, numgp, radius, move, move_by_xx, reflect, importing):
        self.planar     = planar
        self.sphere     = not(self.planar)
        self.grafting   = grafting
        self.radius     = radius
        self.move       = move
        self.move_by_xx = move_by_xx
        self.reflect    = reflect
        self.importing  = importing

        if (self.planar):
            print("Geometry type = planar")
        elif (self.sphere):
            print("Geometry type = spherical")

        if (self.grafting):
            self.numgp       = numgp
            self.graftcoords = np.zeros((3,self.numgp),float)
            self.graftpoints = GraftPoints(self.planar, self.numgp, self.radius)

            if (self.importing):
                self.graftcoords = self.graftpoints.importt()
            else:
                self.graftcoords = self.graftpoints.generate(self.radius)

            self.numgp    = self.graftpoints.numgp
            self.geometry = Geometry(self.graftcoords[0,:], self.graftcoords[1,:], self.graftcoords[2,:], self.move_by_xx)

            if (self.move):
                self.graftcoords[0,:] = self.geometry.moveAlongX()

            if (self.reflect):
                (tempX, tempY, tempZ) = self.geometry.reflectX()

                self.graftcoords = np.zeros((3,tempX.size))

                self.graftcoords[0,:] = tempX
                self.graftcoords[1,:] = tempY
                self.graftcoords[2,:] = tempZ

                self.numgp = tempX.size


    def read_the_mesh(self):
        try:
            mesh_file = open("mesh.in.mphtxt", 'r')
        except:
            print("ERROR OPENING INPUT FILE mesh.in.mphtxt")
            exit()

        try:
            meshpoints_file = open("tmp_meshpoints.txt", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE tmp_meshpoints.txt")
            exit()

        try:
            elem_con_file = open("tmp_elemcon.txt", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE tmp_elemcon.txt")
            exit()

        for line in mesh_file:
            if ("sdim" in line):
                self.ndm = int(line.split()[0])
                print("Number of dimensions        =", self.ndm)

            if ("number of mesh points" in line):
                self.numnp = int(line.split()[0])
                print("Number of nodes             =", self.numnp)
                self.meshpoint = np.zeros((self.ndm,self.numnp))

            if ("Mesh point coordinates" in line):
                for jj in range(0,self.numnp):
                    coords = mesh_file.readline()
                    meshpoints_file.write(coords)
                    #meshpoints_file.write(str(jj+1) + " " + coords)
                    coords = coords.split()
                    self.meshpoint[0,jj] = float(coords[0])
                    self.meshpoint[1,jj] = float(coords[1])
                    self.meshpoint[2,jj] = float(coords[2])

                self.x_min = np.amin(self.meshpoint,axis=1)[0]
                self.y_min = np.amin(self.meshpoint,axis=1)[1]
                self.z_min = np.amin(self.meshpoint,axis=1)[2]

                self.x_max = np.amax(self.meshpoint,axis=1)[0]
                self.y_max = np.amax(self.meshpoint,axis=1)[1]
                self.z_max = np.amax(self.meshpoint,axis=1)[2]

                self.boxX = self.x_max - self.x_min
                self.boxY = self.y_max - self.y_min
                self.boxZ = self.z_max - self.z_min

            if ("4 # number of nodes per element" in line):
                self.nen = int(line.split()[0])
                print("Number of nodes per element =", self.nen)

                self.numel = int(mesh_file.readline().split()[0])
                print("Number of elements          =", self.numel)

                mesh_file.readline()
                self.global_node_id = np.zeros((self.nen,self.numel))

                for kk in range(0,self.numel):
                    global_node_index = mesh_file.readline()
                    elem_con_file.write(global_node_index)
                    #elem_con_file.write(str(kk+1) + " " + global_node_index)
                    global_node_index = global_node_index.split()

                    for pp in range(0,self.nen):
                        self.global_node_id[pp,kk] = global_node_index[pp]

        mesh_file.close()
        meshpoints_file.close()
        elem_con_file.close()

        return


    def matlab(self, modelName):
        try:
            mFile = open(modelName + ".m", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE " + modelName + ".m")
            exit()

        mFile.write("function " + modelName + '\n')
        mFile.write('\n')

        mFile.write("import com.comsol.model.*\n")
        mFile.write("import com.comsol.model.util.*\n")
        mFile.write('\n')

        mFile.write("model = ModelUtil.create('Model');\n")
        mFile.write('\n')

        mFile.write("model.component.create('comp1', false);\n")
        mFile.write('\n')

        mFile.write("model.component('comp1').geom.create('geom1', 3);\n")
        mFile.write('\n')

        for point in range(0, self.numgp):
            mFile.write("model.component('comp1').geom('geom1').create('p" + str(point+1) + "'," + " 'Point');\n")
            mFile.write("model.component('comp1').geom('geom1').feature('p" + str(point+1) + "').set('p',[" + str(self.graftcoords[0,point]) + "," + str(self.graftcoords[1,point]) + "," + str(self.graftcoords[2,point]) + "]);\n")

        mFile.write("model.component('comp1').geom('geom1').run;\n")
        mFile.write('\n')

        mFile.write("DIR = '/home/cjrevelas/';\n")
        mFile.write("FF  = '" + modelName + "';\n")
        mFile.write('\n')

        mFile.write("filename = strcat(FF,'.mph');\n")
        mFile.write("file     = strcat(DIR,filename);\n")
        mFile.write("model.save(file);\n")

        mFile.close()

        return


    def coords_to_nodes(self):
        try:
            nodeFile = open("in.gnodes", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE in.gnodes")
            exit()

        nodeFile.write("ITEM: TIMESTEP\n")
        nodeFile.write("0\n")
        nodeFile.write("ITEM: NUMBER OF ATOMS\n")
        nodeFile.write(str(self.numgp)+'\n')
        nodeFile.write("ITEM: BOX BOUNDS pp pp pp\n")
        nodeFile.write("%.15f %.15f\n" % (-0.5*self.boxX, 0.5*self.boxX) )
        nodeFile.write("%.15f %.15f\n" % (-0.5*self.boxY, 0.5*self.boxY) )
        nodeFile.write("%.15f %.15f\n" % (-0.5*self.boxZ, 0.5*self.boxZ) )
        nodeFile.write("ITEM: ATOMS id type xu yu zu\n")

        tol     = 1.0e-04
        initial = 1.0e+06

        nodeNumber = 0
        for jj in range(0,self.numnp):

            nodeNumber += 1

            for pp in range(0, self.numgp):
                if np.absolute(self.meshpoint[0,jj]-self.graftcoords[0,pp])<tol and np.absolute(self.meshpoint[1,jj]-self.graftcoords[1,pp])<tol and np.absolute(self.meshpoint[2,jj]-self.graftcoords[2,pp])<tol:
                    nodeFile.write(str(nodeNumber)+ "  %.1f" %(initial) + "  %.15f  %.15f  %.15f " %(self.graftcoords[0,pp], self.graftcoords[1,pp], self.graftcoords[2,pp]) + '\n')
                    break

        nodeFile.close()

        return
#----------------------------------------------------------------------------------------------------------------------------------#
class GraftPoints:
    def __init__(self, planar, numgp, radius):
        self.planar = planar
        self.sphere = not(self.planar)
        self.numgp  = numgp
        self.count  = 0


    def importt(self):
        self.coords = np.array([[],[],[]],float)

        try:
            graftFile = open("in.gpoints", 'r')
        except:
            print("ERROR OPENING INPUT FILE in.gpoints")
            exit()

        for line in graftFile:
            x = float(line.split()[0])
            y = float(line.split()[1])
            z = float(line.split()[2])
            self.coords = np.append(self.coords, [[x], [y], [z]], axis=1)
            self.count += 1

        self.numgp = self.count

        return self.coords


    def generate(self, radius):
        if (self.planar):
            self.coords = np.zeros((3,self.numgp),float)

            periodicity_along_x = True
            periodicity_along_y = True
            periodicity_along_z = False

            dz = 0.2

            gp_dist = np.sqrt(1/self.sigma)

            if (periodicity_along_x):
                numgp_x = int(self.boxX / gp_dist)

                for ii in range(0, numgp_x):
                    for jj in range(numgp_x*ii, numgp_x*(ii+1)):
                        self.coords[0,jj] = float(ii+1) * gp_dist - self.boxX * int( (float(ii+1)*gp_dist) / self.boxX) - gp_dist

            if (periodicity_along_y):
                numgp_y = int(self.boxY / gp_dist)

                for ii in range(0, numgp_x*numgp_y):
                    self.coords[1,ii] = float(ii+1) * gp_dist - self.boxX * int( (float(ii+1)*gp_dist) / self.boxY) - gp_dist

            if (periodicity_along_z):
                numgp_z = int(self.boxZ / gp_dist)
            else:
                self.coords[2,:] = self.boxZ - dz

            return self.coords

        if (self.sphere):
            self.radius = radius

            self.coords = np.array([[],[],[]],float)

            #alpha = (4 * np.pi * self.radius**2)/(self.radius**2 * self.num_chains)
            alpha = 4 * np.pi / self.numgp

            d = np.sqrt(alpha)

            m_theta = int(round(np.pi/d))

            d_theta = np.pi/m_theta
            d_phi   = alpha/d_theta

            print("d_theta/d_phi = %.3f" %(d_theta/d_phi))

            for ii in range(0, m_theta):
                theta = np.pi * (ii + 0.5)/m_theta
                m_phi = int(round(2 * np.pi * np.sin(theta)/d_phi))

                for jj in range(0, m_phi):
                    phi = (2 * np.pi * jj)/m_phi

                    x = self.radius * np.sin(theta) * np.cos(phi)
                    y = self.radius * np.sin(theta) * np.sin(phi)
                    z = self.radius * np.cos(theta)

                    self.coords = np.append(self.coords, [[x], [y], [z]], axis=1)

                    self.count += 1

            self.numgp = self.count

            print("Number of grafted chains =", self.count)
            print("Surface grafting density = %.3f chains/Angstrom2" %(float(self.count)/(4*np.pi*self.radius**2)))
            print("                         = %.2f chains/nm2      " %(float(self.count)/(4*np.pi*self.radius**2)*100))

            return self.coords


    def insert_to_model_file(self):
        print("Copying grafting points to matlab model file.")

        try:
            modelFile = open("rsl3d_mesh_gen.m", 'r+')
        except:
            print("ERROR OPENING MODEL FILE rsl3d_mesh_gen.m")
            exit()

        try:
            graftpointsFile = open("tmp_graftpoints.m", 'r')
        except:
            print("ERROR OPENING GRAFTPOINS FILE tmp_graftpoints.m")
            exit()

        lines = modelFile.readlines()
        for num, line_1 in enumerate(lines):
            if ("INSERT GRAFTING POINTS HERE" in line_1):
                break

        for line_2 in reversed(list(graftpointsFile)):
            if ("Point" in line_2 or "set" in line_2):
                lines.insert(num+1, line_2)

        modelFile.truncate(0)
        modelFile.seek(0)
        modelFile.writelines(lines)

        return
#----------------------------------------------------------------------------------------------------------------------------------#
class Geometry:
    def __init__(self, x, y, z, move_by_xx):
        self.x          = x
        self.y          = y
        self.z          = z
        self.move_by_xx = move_by_xx


    def cutInHalf(self):
        positiveX = np.where(self.x > 0.0)
        self.negX = np.delete(self.x, positiveX)
        self.negY = np.delete(self.y, positiveX)
        self.negZ = np.delete(self.z, positiveX)

        negativeX = np.where(self.x < 0.0)
        self.posX = np.delete(self.x, negativeX)
        self.posY = np.delete(self.y, negativeX)
        self.posZ = np.delete(self.z, negativeX)

        return


    def moveAlongX(self):
        self.x += self.move_by_xx

        return self.x


    def reflectX(self):
        self.x_refl = -self.x

        self.x = np.append(self.x, self.x_refl)
        self.y = np.append(self.y, self.y)
        self.z = np.append(self.z, self.z)

        return (self.x, self.y, self.z)
#----------------------------------------------------------------------------------------------------------------------------------#
def set_variable_value(variable, value, fileName):
    strval = "{:.2f}".format(value)
    cmd = "sed -i 's/" + variable + '/' + str(strval) + "/g' " + fileName
    os.system(cmd)

    return


def run_matlab_model():
    print("Running rsl3d_mesh_gen..")
    os.system("./matlab_run.sh")

    return


def edit_model_size_parameters(boxLx, boxLy, boxLz, r_np_eff, centers):
    print("Editing model size parameters..")
    print("boxLx:   ", boxLx)
    set_variable_value("BOXLX", boxLx, "rsl3d_mesh_gen.m")
    print("boxLy:   ", boxLy)
    set_variable_value("BOXLY", boxLy, "rsl3d_mesh_gen.m")
    print("boxLz:   ", boxLz)
    set_variable_value("BOXLZ", boxLz, "rsl3d_mesh_gen.m")
    print("radius: ", r_np_eff)
    set_variable_value("RNP_EFF", r_np_eff, "rsl3d_mesh_gen.m")

    for ii in range(0,centers.shape[1]):
        print("np:       ", ii+1)
        print("x_center: ", centers[0][ii])
        set_variable_value("X_CENTER_SPH" + str(ii+1), centers[0][ii], "rsl3d_mesh_gen.m")
        print("y_center: ", centers[1][ii])
        set_variable_value("Y_CENTER_SPH" + str(ii+1), centers[1][ii], "rsl3d_mesh_gen.m")
        print("z_center: ", centers[2][ii])
        set_variable_value("Z_CENTER_SPH" + str(ii+1), centers[2][ii], "rsl3d_mesh_gen.m")

    return
#----------------------------------------------------------------------------------------------------------------------------------#
planar = False

Lx = 380.0
Ly = 220.0
Lz = 220.0

num_part = 2
r_np_eff = 24.0
r_gp     = 0.4

h_ss_HS  = 132.8  # if num_part=1, this number does not make a difference

use_gr  = True
importt = True
numgp   = 30     # If importt=true (nonuniform), this number does not make a difference

if num_part==1:
    h_cc = 0
elif num_part==2:
    h_cc = 2*r_np_eff + h_ss_HS

centers       = np.empty((3,num_part), float)
centers[0][0] = h_cc / 2.0
centers[1][0] = 0.0
centers[2][0] = 0.0

if num_part == 1:
    os.system("cp ./template_input_files/rsl3d_mesh_gen_template_sph1.m rsl3d_mesh_gen.m")

    move       = False
    move_by_xx = 0.0
    reflect    = False
elif num_part == 2:
    os.system("cp ./template_input_files/rsl3d_mesh_gen_template_sph2.m rsl3d_mesh_gen.m")

    if (importt):
        move       = False
        move_by_xx = 0
        reflect    = False
    else:
        move       = True
        move_by_xx = h_cc / 2.0
        reflect    = True

    centers[0][1] = -centers[0][0]
    centers[1][1] =  centers[1][0]
    centers[2][1] =  centers[2][0]

mesh = Mesh(planar, use_gr, numgp, r_np_eff + r_gp, move, move_by_xx, reflect, importt)

mesh.matlab("tmp_graftpoints")

mesh.graftpoints.insert_to_model_file()

edit_model_size_parameters(Lx, Ly, Lz, r_np_eff, centers)
#exit()
run_matlab_model()

mesh.read_the_mesh()

mesh.coords_to_nodes()

#os.system("rm rsl3d_mesh_gen.m")

exit()
#----------------------------------------------------------------------------------------------------------------------------------#
