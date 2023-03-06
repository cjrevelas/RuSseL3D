#----------------------------------------------------------------------------------------------------------------------------------#
import numpy as np
import os
import sys
#----------------------------------------------------------------------------------------------------------------------------------#
class Mesh:
    def __init__(self, planar, Lx, Ly, Lz, grafting, numgp, radius, moveXX, move_by_xx, moveYY, move_by_yy, moveZZ, move_by_zz, reflectXX, reflectYY, reflectZZ, importing):
        if (planar):
            print("Geometry type = planar")
        else:
            print("Geometry type = spherical")

        if (grafting):
            self.graftpoints = GraftPoints(numgp)

            if (importing):
                self.graftpoints.coords = self.graftpoints.importt()
            else:
                self.graftpoints.coords = self.graftpoints.generate(planar, Lx, Ly, Lz, radius)

            self.geometry = Geometry(self.graftpoints.coords[0,:], self.graftpoints.coords[1,:], self.graftpoints.coords[2,:])

            if (moveXX):
                self.graftpoints.coords[0,:] = self.geometry.moveAlongX(move_by_xx)

            if (moveYY):
                self.graftpoints.coords[1,:] = self.geometry.moveAlongY(move_by_yy)

            if (moveZZ):
                self.graftpoints.coords[2,:] = self.geometry.moveAlongZ(move_by_zz)

            if (reflectXX):
                (tempX, tempY, tempZ) = self.geometry.reflectX()

                self.graftpoints.coords = np.zeros((3,tempX.size))
                self.graftpoints.numgp  = tempX.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ

            if (reflectYY):
                (tempX, tempY, tempZ) = self.geometry.reflectY()

                self.graftpoints.coords = np.zeros((3,tempY.size))
                self.graftpoints.numgp  = tempY.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ

            if (reflectZZ):
                (tempX, tempY, tempZ) = self.geometry.reflectZ()

                self.graftpoints.coords = np.zeros((3,tempZ.size))
                self.graftpoints.numgp  = tempZ.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ


    def read_the_mesh(self, Lx, Ly, Lz):
        try:
            mesh_file = open("in.mesh.mphtxt", 'r')
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
                    coords = coords.split()
                    self.meshpoint[0,jj] = float(coords[0])
                    self.meshpoint[1,jj] = float(coords[1])
                    self.meshpoint[2,jj] = float(coords[2])

                x_min = np.amin(self.meshpoint,axis=1)[0]
                y_min = np.amin(self.meshpoint,axis=1)[1]
                z_min = np.amin(self.meshpoint,axis=1)[2]

                x_max = np.amax(self.meshpoint,axis=1)[0]
                y_max = np.amax(self.meshpoint,axis=1)[1]
                z_max = np.amax(self.meshpoint,axis=1)[2]

                self.boxX = x_max - x_min
                self.boxY = y_max - y_min
                self.boxZ = z_max - z_min

                if (np.absolute(Lx-self.boxX)>1.e-6 or np.absolute(Ly-self.boxY)>1.e-6 or np.absolute(Lz-self.boxZ)>1.e-6):
                    message = "Box dimensions do not match:\n boxX = %.3f, boxY = %.3f, boxZ = %.3f\n Lx = %.3f, Ly = %.3f, Lz = %.3f" % (self.boxX, self.boxY, self.boxZ, Lx, Ly, Lz)
                    raise Exception(message)

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

        for point in range(0, self.graftpoints.numgp):
            mFile.write("model.component('comp1').geom('geom1').create('p" + str(point+1) + "'," + " 'Point');\n")
            mFile.write("model.component('comp1').geom('geom1').feature('p" + str(point+1) + "').set('p',[" + str(self.graftpoints.coords[0,point]) +
                                                                                                        "," + str(self.graftpoints.coords[1,point]) +
                                                                                                        "," + str(self.graftpoints.coords[2,point]) + "]);\n")

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
        nodeFile.write(str(self.graftpoints.numgp)+'\n')
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

            for pp in range(0, self.graftpoints.numgp):
                if (np.absolute(self.meshpoint[0,jj]-self.graftpoints.coords[0,pp])<tol and
                    np.absolute(self.meshpoint[1,jj]-self.graftpoints.coords[1,pp])<tol and
                    np.absolute(self.meshpoint[2,jj]-self.graftpoints.coords[2,pp])<tol):
                    nodeFile.write(str(nodeNumber)+ "  %.1f" %(initial) + "  %.15f  %.15f  %.15f " %(self.graftpoints.coords[0,pp],
                                                                                                     self.graftpoints.coords[1,pp],
                                                                                                     self.graftpoints.coords[2,pp]) + '\n')
                    break

        nodeFile.close()

        return
#----------------------------------------------------------------------------------------------------------------------------------#
class GraftPoints:
    def __init__(self, numgp):
        self.numgp  = numgp
        self.coords = np.zeros((3,self.numgp),float)


    def importt(self):
        try:
            graftFile = open("in.gpoints", 'r')
        except:
            print("ERROR OPENING INPUT FILE in.gpoints")
            exit()

        count = 0
        for line in graftFile:
            self.coords[0,count] = float(line.split()[0])
            self.coords[1,count] = float(line.split()[1])
            self.coords[2,count] = float(line.split()[2])
            count += 1

        self.numgp = count

        return self.coords


    def generate(self, planar, Lx, Ly, Lz, radius):
        if (planar):
            periodicity_along_x = True
            periodicity_along_y = True
            periodicity_along_z = False

            dz = 0.2
            sigma = self.numgp / (Lx * Ly)
            gp_dist = np.sqrt(1/sigma)

            if (periodicity_along_x):
                numgp_x = int(Lx / gp_dist)

                for ii in range(0, numgp_x):
                    for jj in range(numgp_x*ii, numgp_x*(ii+1)):
                        self.coords[0,jj] = float(ii+1) * gp_dist - Lx * int( (float(ii+1)*gp_dist) / Lx) - gp_dist

            if (periodicity_along_y):
                numgp_y = int(Ly / gp_dist)

                for ii in range(0, numgp_x*numgp_y):
                    self.coords[1,ii] = float(ii+1) * gp_dist - Lx * int( (float(ii+1)*gp_dist) / Ly) - gp_dist

            if (not periodicity_along_z):
                self.coords[2,:] = Lz - dz

        if (not(planar)):
            #alpha = (4 * np.pi * radius**2)/(radius**2 * self.numgp)
            alpha = 4 * np.pi / self.numgp

            d = np.sqrt(alpha)

            m_theta = int(round(np.pi/d))

            d_theta = np.pi/m_theta
            d_phi   = alpha/d_theta

            print("d_theta/d_phi = %.3f" %(d_theta/d_phi))

            count = 0
            for ii in range(0, m_theta):
                theta = np.pi * (ii + 0.5)/m_theta
                m_phi = int(round(2 * np.pi * np.sin(theta)/d_phi))

                for jj in range(0, m_phi):
                    phi = (2 * np.pi * jj)/m_phi

                    if (count < self.numgp):
                        self.coords[0,count] = radius * np.sin(theta) * np.cos(phi)
                        self.coords[1,count] = radius * np.sin(theta) * np.sin(phi)
                        self.coords[2,count] = radius * np.cos(theta)
                    else:
                        x = radius * np.sin(theta) * np.cos(phi)
                        y = radius * np.sin(theta) * np.sin(phi)
                        z = radius * np.cos(theta)
                        tempArray = np.array([[x],[y],[z]])
                        self.coords = np.append(self.coords, tempArray, axis=1)

                    count += 1

            self.numgp = count

            print("Number of grafted chains =", self.numgp)
            print("Surface grafting density = %.3f chains/Angstrom2" %(float(self.numgp)/(4*np.pi*radius**2)))
            print("                         = %.2f chains/nm2      " %(float(self.numgp)/(4*np.pi*radius**2)*100))

        return self.coords


    def checkOutOfBox(self, Lx, Ly, Lz):
        self.coords = np.delete(self.coords, np.where(self.coords>=Lx/2.0)[1], axis=1)
        self.coords = np.delete(self.coords, np.where(self.coords<=-Lx/2.0)[1], axis=1)
        self.numgp = self.coords.shape[1]
        for ii in range(self.numgp):
            sys.stdout.write("%d  %.4f  %.4f  %.4f\n" %(ii+1, self.coords[0,ii], self.coords[1,ii], self.coords[2,ii]))
        return


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


    def exportt(self, Lx, Ly, Lz):
        try:
            testFile = open("test.lammpstrj", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE test.lammpstrj")
            exit()

        testFile.write("ITEM: TIMESTEP\n")
        testFile.write("0\n")
        testFile.write("ITEM: NUMBER OF ATOMS\n")
        testFile.write(str(self.numgp)+'\n')
        testFile.write("ITEM: BOX BOUNDS pp pp pp\n")
        testFile.write("%.15f %.15f\n" % (-0.5*Lx, 0.5*Lx) )
        testFile.write("%.15f %.15f\n" % (-0.5*Ly, 0.5*Ly) )
        testFile.write("%.15f %.15f\n" % (-0.5*Lz, 0.5*Lz) )
        testFile.write("ITEM: ATOMS id type xu yu zu\n")

        atomType = 0

        nodeNumber = 0
        for jj in range(0,self.numgp):
            nodeNumber += 1
            testFile.write(str(nodeNumber)+ "  %d" %(atomType) + "  %.15f  %.15f  %.15f " %(self.coords[0,jj],
                                                                                            self.coords[1,jj],
                                                                                            self.coords[2,jj]) + '\n')

        testFile.close()

        return
#----------------------------------------------------------------------------------------------------------------------------------#
class Geometry:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


    def cutInHalf(self):
        positiveX = np.where(self.x > 0.0)
        negX = np.delete(self.x, positiveX)
        negY = np.delete(self.y, positiveX)
        negZ = np.delete(self.z, positiveX)

        negativeX = np.where(self.x < 0.0)
        posX = np.delete(self.x, negativeX)
        posY = np.delete(self.y, negativeX)
        posZ = np.delete(self.z, negativeX)

        return (posX, posY, posZ, negX, negY, negZ)


    def moveAlongX(self, move_by_xx):
        self.x += move_by_xx

        return self.x


    def moveAlongY(self, move_by_yy):
        self.y += move_by_yy

        return self.y


    def moveAlongZ(self, move_by_zz):
        self.z += move_by_zz

        return self.z


    def reflectX(self):
        x_refl = -self.x

        self.x = np.append(self.x, x_refl)
        self.y = np.append(self.y, self.y)
        self.z = np.append(self.z, self.z)

        return (self.x, self.y, self.z)


    def reflectY(self):
        y_refl = -self.y

        self.x = np.append(self.x, self.x)
        self.y = np.append(self.y, y_refl)
        self.z = np.append(self.z, self.z)

        return (self.x, self.y, self.z)


    def reflectZ(self):
        z_refl = -self.z

        self.x = np.append(self.x, self.x)
        self.y = np.append(self.y, self.y)
        self.z = np.append(self.z, z_refl)

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


def edit_model_size_parameters(Lx, Ly, Lz, r_np_eff, centers):
    print("Editing model size parameters..")
    print("boxLx:   ", Lx)
    set_variable_value("BOXLX", Lx, "rsl3d_mesh_gen.m")
    print("boxLy:   ", Ly)
    set_variable_value("BOXLY", Ly, "rsl3d_mesh_gen.m")
    print("boxLz:   ", Lz)
    set_variable_value("BOXLZ", Lz, "rsl3d_mesh_gen.m")
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

Lx = 80.0   # 380.0
Ly = 80.0   # 220.0
Lz = 80.0   # 220.0

num_part = 8  # 2
r_np_eff = 24.0
r_gp     = 0.4

use_gr  = True
importt = False #True
numgp   = 37    # 30 # If importt=true (e.g., nonuniform), this number does not make a difference

centers = np.zeros((3,num_part), float)

if (num_part == 1):
    os.system("cp ./template_input_files/rsl3d_mesh_gen_template_sph1.m rsl3d_mesh_gen.m")

    moveXX     = False
    move_by_xx = 0.0
    reflect    = False
elif (num_part == 2):
    os.system("cp ./template_input_files/rsl3d_mesh_gen_template_sph2.m rsl3d_mesh_gen.m")

    h_ss_HS = 132.8
    h_cc    = 2*r_np_eff + h_ss_HS

    centers[0][0] = +h_cc / 2.0
    centers[0][1] = -h_cc / 2.0

    if (importt):
        moveXX     = False
        move_by_xx = 0
        reflect    = False
    else:
        moveXX     = True
        move_by_xx = h_cc / 2.0
        reflect    = True
else:
    os.system("cp ./template_input_files/rsl3d_mesh_gen_template_sph8.m rsl3d_mesh_gen.m")

    moveXX     = True
    moveYY     = True
    moveZZ     = True
    move_by_xx = Lx / 2.0
    move_by_yy = Ly / 2.0
    move_by_zz = Lz / 2.0
    reflectXX  = True
    reflectYY  = True
    reflectZZ  = True

    centers[0,0] = -Lx / 2.0
    centers[1,0] = -Ly / 2.0
    centers[2,0] = -Lz / 2.0

    centers[0,1] = +Lx / 2.0
    centers[1,1] = -Ly / 2.0
    centers[2,1] = -Lz / 2.0

    centers[0,2] = -Lx / 2.0
    centers[1,2] = +Ly / 2.0
    centers[2,2] = -Lz / 2.0

    centers[0,3] = -Lx / 2.0
    centers[1,3] = -Ly / 2.0
    centers[2,3] = +Lz / 2.0

    centers[0,4] = -Lx / 2.0
    centers[1,4] = +Ly / 2.0
    centers[2,4] = +Lz / 2.0

    centers[0,5] = +Lx / 2.0
    centers[1,5] = -Ly / 2.0
    centers[2,5] = +Lz / 2.0

    centers[0,6] = +Lx / 2.0
    centers[1,6] = +Ly / 2.0
    centers[2,6] = -Lz / 2.0

    centers[0,7] = +Lx / 2.0
    centers[1,7] = +Ly / 2.0
    centers[2,7] = +Lz / 2.0

mesh = Mesh(planar,
            Lx,
            Ly,
            Lz,
            use_gr,
            numgp,
            r_np_eff + r_gp,
            moveXX,
            move_by_xx,
            moveYY,
            move_by_yy,
            moveZZ,
            move_by_zz,
            reflectXX,
            reflectYY,
            reflectZZ,
            importt)

mesh.graftpoints.checkOutOfBox(Lx, Ly, Lz)
#mesh.matlab("tmp_graftpoints")
#mesh.graftpoints.exportt(Lx, Ly, Lz)
#mesh.graftpoints.insert_to_model_file()

#edit_model_size_parameters(Lx, Ly, Lz, r_np_eff, centers)
#exit()
run_matlab_model()

mesh.read_the_mesh(Lx, Ly, Lz)

mesh.coords_to_nodes()

#os.system("rm rsl3d_mesh_gen.m")

exit()
#----------------------------------------------------------------------------------------------------------------------------------#
