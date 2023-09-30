#----------------------------------------------------------------------------------------------------------------------------------#
import numpy as np
import os
import sys
#----------------------------------------------------------------------------------------------------------------------------------#
class Mesh:
    def __init__(self, planar, Lx, Ly, Lz, grafting, numberOfPoints, radius, move, moveBy, reflect, importing, generalPeriodicMesh):

        if (grafting):
            self.graftpoints = GraftPoints(numberOfPoints)

            if (importing):
                self.graftpoints.coords = self.graftpoints.Import(generalPeriodicMesh)
            else:
                self.graftpoints.coords = self.graftpoints.Generate(planar, Lx, Ly, Lz, radius)

            self.geometry = Geometry(self.graftpoints.coords[0,:], self.graftpoints.coords[1,:], self.graftpoints.coords[2,:])

            if (move[0]):
                # TODO: Remove movealong function completely
                self.graftpoints.coords[0,:] = self.geometry.MoveALongX(moveBy[0])

            if (move[1]):
                self.graftpoints.coords[1,:] = self.geometry.MoveALongY(moveBy[1])

            if (move[2]):
                self.graftpoints.coords[2,:] = self.geometry.MoveALongZ(moveBy[2])

            if (reflect[0]):
                (tempX, tempY, tempZ) = self.geometry.ReflectX()

                self.graftpoints.coords = np.zeros((3,tempX.size))
                self.graftpoints.numberOfPoints  = tempX.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ

            if (reflect[1]):
                (tempX, tempY, tempZ) = self.geometry.ReflectY()

                self.graftpoints.coords = np.zeros((3,tempY.size))
                self.graftpoints.numberOfPoints  = tempY.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ

            if (reflect[2]):
                (tempX, tempY, tempZ) = self.geometry.ReflectZ()

                self.graftpoints.coords = np.zeros((3,tempZ.size))
                self.graftpoints.numberOfPoints  = tempZ.size

                self.graftpoints.coords[0,:] = tempX
                self.graftpoints.coords[1,:] = tempY
                self.graftpoints.coords[2,:] = tempZ


    def ReadMesh(self, Lx, Ly, Lz):
        try:
            mesh_file = open("in.mesh", 'r')
        except:
            print("ERROR OPENING INPUT FILE in.mesh")
            exit()

        try:
            meshpoints_file = open("t.meshpoints", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE t.meshpoints")
            exit()

        try:
            elem_con_file = open("t.elemcon", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE t.elemcon")
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


    def Matlab(self, modelName):
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

        for point in range(0, self.graftpoints.numberOfPoints):
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


    def CoordsToNodes(self):
        try:
            nodeFile = open("in.gnodes", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE in.gnodes")
            exit()

        nodeFile.write("ITEM: TIMESTEP\n")
        nodeFile.write("0\n")
        nodeFile.write("ITEM: NUMBER OF ATOMS\n")
        nodeFile.write(str(self.graftpoints.numberOfPoints)+'\n')
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

            for pp in range(0, self.graftpoints.numberOfPoints):
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
    def __init__(self, numberOfPoints):
        self.numberOfPoints  = numberOfPoints
        self.coords = np.zeros((3,self.numberOfPoints),float)


    def Import(self, generalPeriodicMesh):
        if generalPeriodicMesh:

            try:
                graftFile = open("o.mesh_input", 'r')
            except:
                print("ERROR OPENING INPUT FILE o.mesh_input")
                exit()

            count = 0
            for line in graftFile:
                if ("N_GP" in line):
                    count = int(graftFile.readline())

                    self.numberOfPoints = count
                    self.coords = np.zeros((3,self.numberOfPoints),float)

                    graftFile.readline()

                    for point in range(count):
                        nextLine = graftFile.readline()
                        self.coords[0,point] = nextLine.split()[0]
                        self.coords[1,point] = nextLine.split()[1]
                        self.coords[2,point] = nextLine.split()[2]

        else:
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

        self.numberOfPoints = count

        return self.coords


    def Generate(self, planar, Lx, Ly, Lz, radius):
        if (planar):
            periodicity_along_x = True
            periodicity_along_y = True
            periodicity_along_z = False

            dz = 0.2
            sigma = self.numberOfPoints / (Lx * Ly)
            gp_dist = np.sqrt(1/sigma)

            if (periodicity_along_x):
                numberOfPoints_x = int(Lx / gp_dist)

                for ii in range(0, numberOfPoints_x):
                    for jj in range(numberOfPoints_x*ii, numberOfPoints_x*(ii+1)):
                        self.coords[0,jj] = float(ii+1) * gp_dist - Lx * int( (float(ii+1)*gp_dist) / Lx) - gp_dist

            if (periodicity_along_y):
                numberOfPoints_y = int(Ly / gp_dist)

                for ii in range(0, numberOfPoints_x*numberOfPoints_y):
                    self.coords[1,ii] = float(ii+1) * gp_dist - Lx * int( (float(ii+1)*gp_dist) / Ly) - gp_dist

            if (not periodicity_along_z):
                self.coords[2,:] = Lz - dz

        if (not(planar)):
            #alpha = (4 * np.pi * radius**2)/(radius**2 * self.numberOfPoints)
            alpha = 4 * np.pi / self.numberOfPoints

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

                    if (count < self.numberOfPoints):
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

            self.numberOfPoints = count

            print("Number of grafted chains =", self.numberOfPoints)
            print("Surface grafting density = %.3f chains/Angstrom2" %(float(self.numberOfPoints)/(4*np.pi*radius**2)))
            print("                         = %.2f chains/nm2      " %(float(self.numberOfPoints)/(4*np.pi*radius**2)*100))

        return self.coords


    def CheckOutOfBox(self, Lx, Ly, Lz):
        self.coords = np.delete(self.coords, np.where(self.coords>=Lx/2.0)[1], axis=1)
        self.coords = np.delete(self.coords, np.where(self.coords<=-Lx/2.0)[1], axis=1)
        self.numberOfPoints = self.coords.shape[1]
        for ii in range(self.numberOfPoints):
            sys.stdout.write("%d  %.4f  %.4f  %.4f\n" %(ii+1, self.coords[0,ii], self.coords[1,ii], self.coords[2,ii]))
        return


    def InsertToModelFile(self):
        print("Copying grafting point coordinates to matlab model file..")

        try:
            modelFile = open("model.m", 'r+')
        except:
            print("ERROR OPENING MODEL FILE model.m")
            exit()

        try:
            graftpointsFile = open("t.graftpoints.m", 'r')
        except:
            print("ERROR OPENING GRAFTPOINS FILE t.graftpoints.m")
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


    def Export(self, Lx, Ly, Lz):
        try:
            testFile = open("test.lammpstrj", 'w')
        except:
            print("ERROR OPENING OUTPUT FILE test.lammpstrj")
            exit()

        testFile.write("ITEM: TIMESTEP\n")
        testFile.write("0\n")
        testFile.write("ITEM: NUMBER OF ATOMS\n")
        testFile.write(str(self.numberOfPoints)+'\n')
        testFile.write("ITEM: BOX BOUNDS pp pp pp\n")
        testFile.write("%.15f %.15f\n" % (-0.5*Lx, 0.5*Lx) )
        testFile.write("%.15f %.15f\n" % (-0.5*Ly, 0.5*Ly) )
        testFile.write("%.15f %.15f\n" % (-0.5*Lz, 0.5*Lz) )
        testFile.write("ITEM: ATOMS id type xu yu zu\n")

        atomType = 0

        nodeNumber = 0
        for jj in range(0,self.numberOfPoints):
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


    def CutInHalf(self):
        positiveX = np.where(self.x > 0.0)
        negX = np.delete(self.x, positiveX)
        negY = np.delete(self.y, positiveX)
        negZ = np.delete(self.z, positiveX)

        negativeX = np.where(self.x < 0.0)
        posX = np.delete(self.x, negativeX)
        posY = np.delete(self.y, negativeX)
        posZ = np.delete(self.z, negativeX)

        return (posX, posY, posZ, negX, negY, negZ)


    def MoveALongX(self, moveByXX):
        self.x += moveByXX

        return self.x


    def MoveALongY(self, moveByYY):
        self.y += moveByYY

        return self.y


    def MoveALongZ(self, moveByZZ):
        self.z += moveByZZ

        return self.z


    def ReflectX(self):
        x_refl = -self.x

        self.x = np.append(self.x, x_refl)
        self.y = np.append(self.y, self.y)
        self.z = np.append(self.z, self.z)

        return (self.x, self.y, self.z)


    def ReflectY(self):
        y_refl = -self.y

        self.x = np.append(self.x, self.x)
        self.y = np.append(self.y, y_refl)
        self.z = np.append(self.z, self.z)

        return (self.x, self.y, self.z)


    def ReflectZ(self):
        z_refl = -self.z

        self.x = np.append(self.x, self.x)
        self.y = np.append(self.y, self.y)
        self.z = np.append(self.z, z_refl)

        return (self.x, self.y, self.z)
#----------------------------------------------------------------------------------------------------------------------------------#
def SetVariableValue(variable, value, fileName):
    strval = "{:.2f}".format(value)
    cmd = "sed -i 's/" + variable + '/' + str(strval) + "/g' " + fileName
    os.system(cmd)

    return


def RunMatlabModel():
    print("Running model..")
    os.system("./matlab_run.sh")

    return


def EditModelSizeParameters(Lx, Ly, Lz, radiusEffective, centers):
    print("Editing model size parameters..")

    SetVariableValue("BOXLX", Lx, "model.m")
    SetVariableValue("BOXLY", Ly, "model.m")
    SetVariableValue("BOXLZ", Lz, "model.m")
    SetVariableValue("RNP_EFF", radiusEffective, "model.m")

    for ii in range(0,centers.shape[1]):
        SetVariableValue("X_CENTER_SPH" + str(ii+1), centers[0][ii], "model.m")
        SetVariableValue("Y_CENTER_SPH" + str(ii+1), centers[1][ii], "model.m")
        SetVariableValue("Z_CENTER_SPH" + str(ii+1), centers[2][ii], "model.m")

    return


def InsertNanoparticlesToModelFile(numberOfParticles, effectiveRadius, denseMeshWidth, centers):
    print("Copying nanoparticle center coordinates to matlab model file..")

    try:
        modelFile = open("model.m", 'r+')
    except:
        print("ERROR OPENING MODEL FILE model.m")
        exit()

    lines = modelFile.readlines()
    for num, line in enumerate(lines):
        if ("INSERT NANOPARTICLES HERE" in line):
            break

    for ii in reversed(range(numberOfParticles)):
        newLine4 = "model.component('comp1').geom('geom1').create('sph"  + str(numberOfParticles+ii+1) + "', 'Sphere');\n"
        newLine5 = "model.component('comp1').geom('geom1').feature('sph" + str(numberOfParticles+ii+1) + "').set('pos',[" + str(centers[0,ii]) + " " + str(centers[1,ii]) + " " + str(centers[2,ii]) + "]);\n"
        newLine6 = "model.component('comp1').geom('geom1').feature('sph" + str(numberOfParticles+ii+1) + "').set('r'," + str(effectiveRadius + denseMeshWidth) + ");\n"

        lines.insert(num+1, newLine6)
        lines.insert(num+1, newLine5)
        lines.insert(num+1, newLine4)

    for ii in reversed(range(numberOfParticles)):
        newLine1 = "model.component('comp1').geom('geom1').create('sph"  + str(ii+1) + "', 'Sphere');\n"
        newLine2 = "model.component('comp1').geom('geom1').feature('sph" + str(ii+1) + "').set('pos',[" + str(centers[0,ii]) + " " + str(centers[1,ii]) + " " + str(centers[2,ii]) + "]);\n"
        newLine3 = "model.component('comp1').geom('geom1').feature('sph" + str(ii+1) + "').set('r'," + str(effectiveRadius) + ");\n"

        lines.insert(num+1, newLine3)
        lines.insert(num+1, newLine2)
        lines.insert(num+1, newLine1)


    modelFile.truncate(0)
    modelFile.seek(0)
    modelFile.writelines(lines)

    return
#----------------------------------------------------------------------------------------------------------------------------------#
# Initialize all variables
generalPeriodicMesh = True

planar     = False
useGrafted = True
importt    = True

numberOfParticles = 14

centers = np.zeros((3,numberOfParticles), float)

Lx = 0.0
Ly = 0.0
Lz = 0.0

radiusEffective = 44.0
numberOfPoints  = 0

graftingDistance = 0.4
denseMeshWidth   = 5.0

move    = np.array([False, False, False], bool)
moveBy  = np.array([False, False, False], bool)
reflect = np.array([False, False, False], bool)

# Check different cases
if (numberOfParticles == 1):
    os.system("cp ./template_input_files/rsl3d_gen_mesh_template_sph1.m model.m")

elif (numberOfParticles == 2):
    os.system("cp ./template_input_files/rsl3d_gen_mesh_template_sph2.m model.m")

    h_ss_HS = 132.8
    h_cc    = 2*radiusEffective + h_ss_HS

    centers[0][0] = + h_cc / 2.0
    centers[0][1] = - h_cc / 2.0

    if (not importt):
        move[0]    = True
        moveBy[0]  = h_cc / 2.0
        reflect[0] = True

elif (generalPeriodicMesh):
    os.system("cp rsl3d_gen_mesh_template_fcc.m model.m")

    try:
        inputFile = open("o.mesh_input", 'r')
    except:
        print("ERROR OPENING INPUT FILE o.mesh_input")
        exit()

    for line in inputFile:
        if ("BOX" in line):
            nextLine = inputFile.readline()
            xlo = float(nextLine.split()[0])
            xhi = float(nextLine.split()[1])
            Lx  = xhi - xlo

            nextLine = inputFile.readline()
            ylo = float(nextLine.split()[0])
            yhi = float(nextLine.split()[1])
            Ly  = yhi - ylo

            nextLine = inputFile.readline()
            zlo = float(nextLine.split()[0])
            zhi = float(nextLine.split()[1])
            Lz  = zhi - zlo

        if ("N_NP" in line):
            numberOfParticles = int(inputFile.readline())

            nextLine = inputFile.readline()

            for particle in range(numberOfParticles):
                nextLine = inputFile.readline()
                centers[0,particle] = float(nextLine.split()[0])
                centers[1,particle] = float(nextLine.split()[1])
                centers[2,particle] = float(nextLine.split()[2])

    InsertNanoparticlesToModelFile(numberOfParticles, radiusEffective, denseMeshWidth, centers)

#Read mesh and generate grafting points
mesh = Mesh(planar, Lx, Ly, Lz, useGrafted, numberOfPoints, radiusEffective + graftingDistance, move, moveBy, reflect, importt, generalPeriodicMesh)

#mesh.graftpoints.CheckOutOfBox(Lx, Ly, Lz)
mesh.Matlab("t.graftpoints")
#mesh.graftpoints.Export(Lx, Ly, Lz)
mesh.graftpoints.InsertToModelFile()

EditModelSizeParameters(Lx, Ly, Lz, radiusEffective, centers)
RunMatlabModel()

print("Reading mesh created by comsol model..")
os.system("mv in.mesh.mphtxt in.mesh")
mesh.ReadMesh(Lx, Ly, Lz)

print("Searching for grafting points for in.gnodes file..")
mesh.CoordsToNodes()

print("Removing temporary files..")
#os.system("rm t.meshpoints t.elemcon")
#os.system("rm model.m")
#os.system("rm t.graftpoints.m")

exit()
#----------------------------------------------------------------------------------------------------------------------------------#
