#!/home/e89/e89/zd242/src/anaconda2/bin/python
# Application for mechanical properties calculation and analysis using VASP
# Zeyu Deng 16.03.2017

import os
import sys
import shutil
import subprocess
import numpy as np
import argparse
import sympy
#import matplotlib.pyplot as plt
###################################################### Constants #########
# Info
info_text = "Mechanical Properties Calculation and Analysis Application ver 1.05\n" +\
    "Zeyu Deng <zd242@cam.ac.uk or dengzeyu@gmail.com>\n" +\
    "Department of Materials Science and Metallurgy\n" +\
    "University of Cambridge\n" +\
    "16.03.2017"


def printInfo():
    print info_text

# POSCAR
poscar = "POSCAR"
# Files to copy to the calculation directories
file_list = [poscar]
c11,c22,c33,c44,c55,c66,c12,c13,c14,c15,c16,c23,c24,c25,c26,c34,c35,c36,c45,c46,c56,d = sympy.symbols("c11,c22,c33,c44,c55,c66,c12,c13,c14,c15,c16,c23,c24,c25,c26,c34,c35,c36,c45,c46,c56,d")

def get_strain_pattern(crystSys):
    if crystSys == 1:  # cubic -> e1+e4
        pattern = [[1, 0, 0, 1, 0, 0]]

    # tetragonal high 4mm, -42m, 422, 4/mmm; tetragonal low 4, -4, 4/m ->
    # e1+e4, e3+e6
    elif crystSys in [2, 21]:
        pattern = [[1, 0, 0, 1, 0, 0],
                   [0, 0, 1, 0, 0, 1]]

    elif crystSys == 3:  # orthohombic -> e1+e4, e2+e5, e3+e6
        pattern = [[1, 0, 0, 1, 0, 0],
                   [0, 1, 0, 0, 1, 0],
                   [0, 0, 1, 0, 0, 1]]

    elif crystSys == 4:  # monoclinic type I Diad || x2 -> e1+e4, e3+e6, e2, e5
        pattern = [[1, 0, 0, 1, 0, 0],
                   [0, 0, 1, 0, 0, 1],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1, 0]]

    # monoclinic type II Diad || x3 -> e1+e4, e3+e5, e2, e6
    elif crystSys == 41:
        pattern = [[1, 0, 0, 1, 0, 0],
                   [0, 0, 1, 0, 1, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 1]]

    elif crystSys == 5:  # triclinic -> e1, e2, e3, e4, e5, e6
        pattern = [[1, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 1]]

    # hexagonal; trignal high 32, -3m, 3m; trignal low 3, -3 -> e1, e3+e4
    elif crystSys in [6, 7, 71]:
        pattern = [[1, 0, 0, 0, 0, 0],
                   [0, 0, 1, 1, 0, 0]]

    return pattern

def get_strain_pattern_sym(crystSys):
	return sympy.Matrix(get_strain_pattern(crystSys))*d

def get_stress_pattern_sym(crystSys): # linear algebra calculation to get symbolic expression of stress
    strain_pattern_sym=get_strain_pattern_sym(crystSys)
    cij_pattern_sym = get_cij_pattern_sym(crystSys)
    stress_pattern_sym = sympy.Matrix([[0, 0, 0, 0, 0, 0]])

    for i in range(strain_pattern_sym.shape[0]):
        stress_pattern_sym = stress_pattern_sym.row_insert(-1, (cij_pattern_sym*strain_pattern_sym.row(i).T).T)

    stress_pattern_sym.row_del(-1)

    return sympy.simplify(stress_pattern_sym)


def get_cij_pattern_sym(crystSys):  # convert cij_pattern into symbolic matrix
    cij_to_string = {1: 'c11', 2: 'c22', 3: 'c33', 4: 'c44', 5: 'c55', 6: 'c66',
                     7: 'c12', 8: 'c13', 9: 'c14', 10: 'c15', 11: 'c16',
                     12: 'c23', 13: 'c24', 14: 'c25', 15: 'c26',
                     16: 'c34', 17: 'c35', 18: 'c36',
                     19: 'c45', 20: 'c46',
                     21: 'c56', 0: 0, -9: '-c14', -14: '-c25', -11: '-c16'}
    cij_pattern = get_cij_pattern(crystSys)
    cij_pattern_sym = []
    for index_i, i in enumerate(cij_pattern):
        cij_pattern_sym_row = []
        for index_j, j in enumerate(i):
            cij_pattern_sym_row.append(cij_to_string[cij_pattern[index_i][index_j]])
        cij_pattern_sym.append(cij_pattern_sym_row)
    cij_pattern_sym = sympy.Matrix(cij_pattern_sym)
    return cij_pattern_sym



def get_deform(crystSys, strain_index=0, deform=0):
    # Cubic volume-conserving for energy method
    if crystSys == 100:
        deformation = [[[1+deform, 0, 0], [0, 1+deform, 0], [0, 0, 1/((1+deform)**2)]],
                       [[1+deform, 0, 0], [0, 1+deform, 0], [0, 0, 1+deform]],
                       [[1, deform/2, 0], [deform/2, 1, 0], [0, 0, 1+(deform**2)/(4-deform**2)]]]
    # Stress-Strain methods
    if crystSys < 100:
        if crystSys != 0:
            strain = get_strain_pattern(crystSys)
            deformation = np.identity(3)+tran(strain[strain_index])*deform
    return deformation


def get_deform_file():
    f = open("deformation.dat", 'r')
    lines = f.readlines()
    f.close()
    strain_list = [[float(num) for num in line.strip().split()]
                   for line in lines]
    deformation = [np.identity(3)+tran(strain_list[strain_index])
                   for strain_index in range(len(strain_list))]
    return deformation


# See J.F. Nye, Physical Properties of Crystals, P140-141
def get_cij_pattern(crystSys):
    if crystSys == 1:  # cubic -> e1+e4
        pattern = [[1, 7, 7, 0, 0, 0],
                   [7, 1, 7, 0, 0, 0],
                   [7, 7, 1, 0, 0, 0],
                   [0, 0, 0, 4, 0, 0],
                   [0, 0, 0, 0, 4, 0],
                   [0, 0, 0, 0, 0, 4]]

    # tetragonal high 4mm, -42m, 422, 4/mmm -> e1+e4, e3+e6
    elif crystSys == 2:
        pattern = [[1, 7, 8, 0, 0, 0],
                   [7, 1, 8, 0, 0, 0],
                   [8, 8, 3, 0, 0, 0],
                   [0, 0, 0, 4, 0, 0],
                   [0, 0, 0, 0, 4, 0],
                   [0, 0, 0, 0, 0, 6]]

    elif crystSys == 21:  # tetragonal low 4, -4, 4/m -> e1+e4, e3+e6
        pattern = [[1, 7, 8, 0, 0, 11],
                   [7, 1, 8, 0, 0, -11],
                   [8, 8, 3, 0, 0, 0],
                   [0, 0, 0, 4, 0, 0],
                   [0, 0, 0, 0, 4, 0],
                   [11, -11, 0, 0, 0, 6]]

    elif crystSys == 3:  # orthohombic -> e1+e4, e2+e5, e3+e6
        pattern = [[1,  7,  8,  0,  0,  0],
                   [7,  2, 12,  0,  0,  0],
                   [8, 12,  3,  0,  0,  0],
                   [0,  0,  0,  4,  0,  0],
                   [0,  0,  0,  0,  5,  0],
                   [0,  0,  0,  0,  0,  6]]

    elif crystSys == 4:  # monoclinic type I Diad || x2 -> e1+e4, e3+e6, e2, e5
        pattern = [[1,  7,  8,  0,  10,  0],
                   [7,  2, 12,  0, 14,  0],
                   [8, 12,  3,  0, 17,  0],
                   [0,  0,  0,  4,  0,  20],
                   [10, 14, 17,  0,  5,  0],
                   [0,  0,  0, 20,  0,  6]]

    # monoclinic type II Diad || x3 -> e1+e4, e3+e5, e2, e6
    elif crystSys == 41:
        pattern = [[1,  7,  8,  0,  0,  11],
                   [7,  2, 12,  0,  0,  15],
                   [8, 12,  3,  0,  0,  18],
                   [0,  0,  0,  4,  19,  0],
                   [0,  0,  0,  19,  5,  0],
                   [11, 15, 18,  0,  0,  6]]

    elif crystSys == 5:  # triclinic -> e1, e2, e3, e4, e5, e6
        pattern = [[1,  7,  8,  9,  10, 11],
                   [7,  2, 12,  13, 14,  15],
                   [8, 12,  3,  16, 17,  18],
                   [9, 13, 16,  4,  19,  20],
                   [10, 14, 17, 19,  5,  21],
                   [11, 15, 18, 20,  21,  6]]

    elif crystSys == 6:  # hexagonal -> e1, e3+e4
        pattern = [[1,  7,  8,  0,  0,  0],
                   [7,  1,  8,  0,  0,  0],
                   [8,  8,  3,  0,  0,  0],
                   [0,  0,  0,  4,  0,  0],
                   [0,  0,  0,  0,  4,  0],
                   [0,  0,  0,  0,  0,  6]]

    elif crystSys == 7:  # trignal high 32, -3m, 3m -> e1, e3+e4
        pattern = [[1,  7,  8,  9,  0,  0],
                   [7,  1,  8, -9,  0,  0],
                   [8,  8,  3,  0,  0,  0],
                   [9, -9,  0,  4,  0,  0],
                   [0,  0,  0,  0,  4,  9],
                   [0,  0,  0,  0,  9,  6]]

    elif crystSys == 71:  # trignal low 3, -3 -> e1, e3+e4
        pattern = [[1,  7,  8,  9, -14,  0],
                   [7,  1,  8, -9,  14,  0],
                   [8,  8,  3,  0,  0,  0],
                   [9, -9,  0,  4,  0, 14],
                   [-14, 14,  0,  0,  4,  9],
                   [0,  0,  0, 14,  9,  6]]
    return pattern


def get_num_strain(crystSys):
    if crystSys == 1:
        return 1
    elif crystSys in [2, 21]:
        return 2
    elif crystSys == 3:
        return 3
    elif crystSys in [4, 41]:
        return 4
    elif crystSys == 5:
        return 6
    elif crystSys in [6, 7, 71]:
        return 2
    return num_strain

voigt_mat = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])

###################################################### Preprocessing Funct


def preprocess(crystSys, num_calcPoint, delta):
    [p, bot, top] = readPoscar()
    deform_matrix = []
    shutil.move(poscar, 'POSCAR_backup')
    num_strain = get_num_strain(crystSys)
    if num_calcPoint % 2 != 0:
        num_calcPoint = num_calcPoint-1
    for index_Strain in range(num_strain):
        deform = -num_calcPoint/2*delta
        mkdir("./str"+str(index_Strain+1))
        for j in range(num_calcPoint/2):
            deform_dir = "./str"+str(index_Strain+1)+"/d"+str(deform)
            mkdir(deform_dir)
            deformation = get_deform(crystSys, index_Strain, deform)
            pos = apply_deform(p, deformation)
            writePoscar(pos, bot, top)
            deform = deform+delta
            copyfiles(file_list, deform_dir)
        deform = delta
        for j in range(num_calcPoint/2):
            deform_dir = "./str"+str(index_Strain+1)+"/d"+str(deform)
            mkdir(deform_dir)
            deformation = get_deform(crystSys, index_Strain, deform)
            pos = apply_deform(p, deformation)
            writePoscar(pos, bot, top)
            deform = deform+delta
            copyfiles(file_list, deform_dir)
        print "Strain complete! (%d/%d)" % (index_Strain+1, num_strain)
        index_Strain = index_Strain+1
    shutil.move('POSCAR_backup', poscar)
    return elastic


def readPoscar():
    f = open(poscar, 'r')
    lines = f.readlines()
    f.close()
    index = 2
    scal = float(lines[index-1].strip().split()[0])
    p = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            p[i][j] = float(lines[3+i-1].strip().split()[j])
    print "POSCAR has been read!"
    print "\n-----------Lattice Vector-------------"
    print scal
    print p
    print "--------------------------------------"
    top = lines[0:2]
    bot = lines[5:len(lines)]
    return p, bot, top


def writePoscar(pos, bot, top, file_name='POSCAR'):
    np.savetxt('poscar', pos, fmt='%21.16f')
    f = open('poscar', 'r')
    lines = f.readlines()
    f.close()
    posstring = top+lines+bot
    pos_file = open(file_name, 'w')
    pos_file.write(''.join(posstring))
    pos_file.close()
    os.remove('poscar')


def copyfiles(file_list, dst):
    for file in file_list:
        try:
            shutil.copy2(file, dst+'/'+file)
        except IOError:
            print "WARNING! Problems on copy file: "+file
            pass


def mkdir(folder):
    try:
        os.mkdir(folder)
    except OSError:
        print "WARNING! File/Folder exists! ("+folder+")"
        pass


def tran(mat):
    mat_tr = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if voigt_mat[i, j] > 2:
                coeff = 0.5
            else:
                coeff = 1
            mat_tr[i][j] = coeff*mat[voigt_mat[i, j]]
    return np.array(mat_tr)


def apply_deform(p, deformation):
    print "\n-------------Deformation--------------"
    print np.array(deformation)
    pos = np.transpose(np.dot(deformation, np.transpose(p)))
    print "----------Deformed Structure----------"
    print np.array(pos)
    print "--------------------------------------"
    return pos

###################################################### Postprocessing Func


def postprocess(arguments):
    crystSys = arguments.crystSys
    c, std_err = fitCij(crystSys)
    Cij = elast_consts(arguments, c, std_err)

def postprocess_test(arguments):
    crystSys = arguments.crystSys
    c, std_err = fitCij(crystSys,arguments)
    Cij = elast_consts(arguments, c, std_err)

def postprocess_read_cij(arguments):
    print "Reading data...."
    c = readCij()
    Cij = elast_consts(arguments, c)


def fitCij(crystSys,arguments):
    print "\n------------------------------------Fitted Results-------------------------------------"

    def fit(i, j):
        from scipy import stats, sqrt, square
        stress_fit = np.array(stress[i, j])
        slope, intercept, r, p, stderr = stats.linregress(delta, stress_fit)
        print '\n'
        print 'Fitted result ', Stress_string[i, j], '          :    ', slope
        print 'Error            	:    ', stderr
        if abs(r) > 0.9:
            print 'Correlation coefficient (r) :    ', r
        else:
            print 'Correlation coefficient (r) :    ', r, '     <----- WARNING: BAD FIT'

        return slope, stderr

    def createCij():
        CijMatrix = np.zeros((6, 6))
        CijErrorMatrix = np.zeros((6, 6))
        c = np.array(get_cij_pattern(crystSys))
        for i in range(0, 6):
            for j in range(0, 6):
                index = int(c[i, j])
                if index > 0:
                    CijMatrix[i, j] = Cij_list[index-1]
                    CijErrorMatrix[i, j] = abs(Cij_errors_list[index-1])
                elif index < 0:
                    CijMatrix[i, j] = -Cij_list[-index-1]
                    CijErrorMatrix[i, j] = abs(Cij_errors_list[-index-1])
        return CijMatrix, CijErrorMatrix

    # get stress from OUTCAR
    command = (
        "grep 'in kB' outcar.2 | tail -1 | awk '{print -$3/10.0, -$4/10.0, -$5/10.0, -$7/10.0, -$8/10.0, -$6/10.0}'")
    stress = []
    strain = []
    delta_list = []
    str_list = os.listdir("./")
    str_list = [int(elem[3:]) for elem in str_list if "str" in elem]
    str_list.sort()
    os.chdir("./rlx")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    stress0 = [float(elem) for elem in process.stdout.read().strip().split()]
    os.chdir("../")
    for this_stress in str_list:
        stress_str = []
        delta_list_in = []
        os.chdir("./str"+str(this_stress))
        delta_string = os.listdir("./")
        delta = [float(elem[1:]) for elem in delta_string if "d" in elem]
        delta.sort()
        for this_delta in delta:
            delta_list_in.append(this_delta)
            os.chdir("./d"+str(this_delta))
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, shell=True)
            stress_str.append(
                [float(elem) for elem in process.stdout.read().strip().split()])
            os.chdir("../")
        os.chdir("../")
        stress.append(stress_str)
        delta_list.append(delta_list_in)
    strain.append(get_strain_pattern(crystSys))
    delta = np.array(delta)

    stress = (np.array(stress)-np.array(stress0)).tolist()
    c = np.zeros((6, 6))
    re_stress = np.zeros((6, 6))
    r_value = np.zeros((6, 6))
    p_value = np.zeros((6, 6))
    std_err = np.zeros((6, 6))
    stress = np.array([np.array(stress_elem).T for stress_elem in stress])

    # Start fitting Cijs
    Stress_string = np.chararray((6, 6), itemsize=15)
    Cij_list = np.zeros(21)
    Cij_errors_list = np.zeros(21)

    print "This code is trying to evaluate: "
    sympy.init_printing()
    stress_sym=(get_stress_pattern_sym(arguments.crystSys)/d)
    strain_sym=get_strain_pattern_sym(arguments.crystSys)
    cij_sym=get_cij_pattern_sym(arguments.crystSys)
    sympy.pprint(sympy.relational.Eq((stress_sym*d).T,sympy.MatMul(cij_sym,strain_sym.T)))

    # Initialize string for Stresses using symbolic notation e.g. C11+C12
    stress_sym=np.array(stress_sym)
    for i in range(stress_sym.shape[0]):
        for j in range(stress_sym.shape[1]):
            Stress_string[i][j]=stress_sym[i][j]

    # Fitting results stored in Fit_results and error in Fit_errors
    Fit_results=np.zeros((stress_sym.shape))
    Fit_errors=np.zeros((stress_sym.shape))

    for i in range(stress_sym.shape[0]):
        for j in range(stress_sym.shape[1]):
            if stress_sym[i][j] != 0:
                Fit_results[i][j],Fit_errors[i][j]= fit(i,j)

    # Solve symbolic equations to get values for each Cijs


    # if crystSys == 1:  # cubic
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C12"
    #     Stress_string[0, 3] = "C44"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit01, error01 = fit(0, 1)  # 12
    #     fit02, error02 = fit(0, 2)  # 12
    #     Cij_list[6] = (fit01+fit02)/2  # c12
    #     Cij_errors_list[6] = (error01+error02)/2

    # elif crystSys == 2:  # tetragonal high 4mm,-42m,422,4/mmm e1+e4, e3+e6
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C44"
    #     Stress_string[1, 0] = "C13"
    #     Stress_string[1, 1] = "C13"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 5] = "C66"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[6], Cij_errors_list[6] = fit(0, 1)  # c12
    #     fit02, error02 = fit(0, 2)  # 13
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit10, error10 = fit(1, 0)  # 13
    #     fit11, error11 = fit(1, 1)  # 13
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     Cij_list[5], Cij_errors_list[5] = fit(1, 5)  # c66
    #     Cij_list[7] = (fit02+fit10+fit11)/3  # c13
    #     Cij_errors_list[7] = (error02+error10+error11)/3

    # elif crystSys == 21:  # tetragonal low 4, -4, 4/m e1+e4, e3+e6
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C44"
    #     Stress_string[0, 5] = "C16"
    #     Stress_string[1, 0] = "C13 + C16"
    #     Stress_string[1, 1] = "C13 - C16"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 5] = "C66"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[6], Cij_errors_list[6] = fit(0, 1)  # c12
    #     fit02, error02 = fit(0, 2)  # 13
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit05, error05 = fit(0, 5)  # 16
    #     fit10, error10 = fit(1, 0)  # 13+16
    #     fit11, error11 = fit(1, 1)  # 13-16
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     Cij_list[5], Cij_errors_list[5] = fit(1, 5)  # c66
    #     Cij_list[7] = (fit02+fit10+fit11)/3  # c13
    #     Cij_errors_list[7] = (error02+error10+error11)/3
    #     Cij_list[10] = (fit10-fit11+fit05)/3  # c16
    #     Cij_errors_list[10] = (error10-error11+error05)/3

    # elif crystSys == 3:  # orthohombic
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C44"
    #     Stress_string[1, 0] = "C12"
    #     Stress_string[1, 1] = "C22"
    #     Stress_string[1, 2] = "C23"
    #     Stress_string[1, 4] = "C55"
    #     Stress_string[2, 0] = "C13"
    #     Stress_string[2, 1] = "C23"
    #     Stress_string[2, 2] = "C33"
    #     Stress_string[2, 5] = "C66"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     fit01, error01 = fit(0, 1)  # 12
    #     fit02, error02 = fit(0, 2)  # 13
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit10, error10 = fit(1, 0)  # 12
    #     Cij_list[1], Cij_errors_list[1] = fit(1, 1)  # c22
    #     fit12, error12 = fit(1, 2)  # 23
    #     Cij_list[4], Cij_errors_list[4] = fit(1, 4)  # c55
    #     fit20, error20 = fit(2, 0)  # 13
    #     fit21, error21 = fit(2, 1)  # 23
    #     Cij_list[2], Cij_errors_list[2] = fit(2, 2)  # c33
    #     Cij_list[5], Cij_errors_list[5] = fit(2, 5)  # c66
    #     Cij_list[6] = (fit10+fit01)/2  # c12
    #     Cij_errors_list[6] = (error10+error01)/2
    #     Cij_list[7] = (fit02+fit20)/2  # c13
    #     Cij_errors_list[7] = (error02+error20)/2
    #     Cij_list[11] = (fit12+fit21)/2  # c23
    #     Cij_errors_list[11] = (error12+error21)/2

    # elif crystSys == 4:  # monoclinic
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C44"
    #     Stress_string[0, 4] = "C15"
    #     Stress_string[0, 5] = "C46"
    #     Stress_string[1, 0] = "C13"
    #     Stress_string[1, 1] = "C23"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 3] = "C46"
    #     Stress_string[1, 4] = "C53"
    #     Stress_string[1, 5] = "C66"
    #     Stress_string[2, 0] = "C12"
    #     Stress_string[2, 1] = "C22"
    #     Stress_string[2, 2] = "C23"
    #     Stress_string[2, 4] = "C25"
    #     Stress_string[3, 0] = "C15"
    #     Stress_string[3, 1] = "C25"
    #     Stress_string[3, 2] = "C35"
    #     Stress_string[3, 4] = "C55"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     fit01, error01 = fit(0, 1)  # 12
    #     fit02, error02 = fit(0, 2)  # 13
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit04, error04 = fit(0, 4)  # 15
    #     fit05, error05 = fit(0, 5)  # 46
    #     fit10, error10 = fit(1, 0)  # 13
    #     fit11, error11 = fit(1, 1)  # 23
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     fit13, error13 = fit(1, 3)  # 46
    #     fit14, error14 = fit(1, 4)  # 35
    #     Cij_list[5], Cij_errors_list[5] = fit(1, 5)  # c66
    #     fit20, error20 = fit(2, 0)  # 12
    #     Cij_list[1], Cij_errors_list[1] = fit(2, 1)  # c22
    #     fit22, error22 = fit(2, 2)  # 23
    #     fit24, error24 = fit(2, 4)  # 25
    #     fit30, error30 = fit(3, 0)  # 15
    #     fit31, error31 = fit(3, 1)  # 25
    #     fit32, error32 = fit(3, 2)  # 35
    #     Cij_list[4], Cij_errors_list[4] = fit(3, 4)  # c55
    #     Cij_list[6] = (fit20+fit01)/2  # c12
    #     Cij_errors_list[6] = (error20+error01)/2
    #     Cij_list[7] = (fit02+fit10)/2  # c13
    #     Cij_errors_list[7] = (error02+error10)/2
    #     Cij_list[11] = (fit11+fit22)/2  # c23
    #     Cij_errors_list[11] = (error11+error22)/2
    #     Cij_list[9] = (fit04+fit30)/2  # c15
    #     Cij_errors_list[9] = (error04+error30)/2
    #     Cij_list[13] = (fit24+fit31)/2  # c25
    #     Cij_errors_list[13] = (error24+error31)/2
    #     Cij_list[16] = (fit14+fit32)/2  # c35
    #     Cij_errors_list[16] = (error14+error32)/2
    #     Cij_list[19] = (fit13+fit05)/2  # c46
    #     Cij_errors_list[19] = (error13+error05)/2

    # elif crystSys == 41:  # monoclinic
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C44"
    #     Stress_string[0, 4] = "C45"
    #     Stress_string[0, 5] = "C16"
    #     Stress_string[1, 0] = "C13"
    #     Stress_string[1, 1] = "C23"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 3] = "C45"
    #     Stress_string[1, 4] = "C55"
    #     Stress_string[1, 5] = "C36"
    #     Stress_string[2, 0] = "C12"
    #     Stress_string[2, 1] = "C22"
    #     Stress_string[2, 2] = "C23"
    #     Stress_string[2, 5] = "C26"
    #     Stress_string[3, 0] = "C16"
    #     Stress_string[3, 1] = "C26"
    #     Stress_string[3, 2] = "C36"
    #     Stress_string[3, 5] = "C66"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     fit01, error01 = fit(0, 1)  # 12
    #     fit02, error02 = fit(0, 2)  # 13
    #     Cij_list[3], Cij_errors_list[3] = fit(0, 3)  # c44
    #     fit04, error04 = fit(0, 4)  # 45
    #     fit05, error05 = fit(0, 5)  # 16
    #     fit10, error10 = fit(1, 0)  # 13
    #     fit11, error11 = fit(1, 1)  # 23
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     fit13, error13 = fit(1, 3)  # 45
    #     Cij_list[4], Cij_errors_list[4] = fit(1, 4)  # c55
    #     fit15, error15 = fit(1, 5)  # 36
    #     fit20, error20 = fit(2, 0)  # 12
    #     Cij_list[1], Cij_errors_list[1] = fit(2, 1)  # c22
    #     fit22, error22 = fit(2, 2)  # 23
    #     fit25, error25 = fit(2, 5)  # 26
    #     fit30, error30 = fit(3, 0)  # 16
    #     fit31, error31 = fit(3, 1)  # 26
    #     fit32, error32 = fit(3, 2)  # 36
    #     Cij_list[5], Cij_errors_list[5] = fit(3, 5)  # c66
    #     Cij_list[6] = (fit20+fit01)/2  # c12
    #     Cij_errors_list[6] = (error20+error01)/2
    #     Cij_list[7] = (fit02+fit10)/2  # c13
    #     Cij_errors_list[7] = (error02+error10)/2
    #     Cij_list[11] = (fit22+fit11)/2  # c23
    #     Cij_errors_list[11] = (error22+error11)/2
    #     Cij_list[10] = (fit05+fit30)/2  # c16
    #     Cij_errors_list[10] = (error05+error30)/2
    #     Cij_list[14] = (fit25+fit31)/2  # c26
    #     Cij_errors_list[14] = (error25+error31)/2
    #     Cij_list[17] = (fit15+fit15)/2  # c36
    #     Cij_errors_list[17] = (error32+error32)/2
    #     Cij_list[18] = (fit13+fit04)/2  # c45
    #     Cij_errors_list[18] = (error13+error04)/2

    # elif crystSys == 5:  # triclinic
    #     index = 0
    #     c = np.array(get_cij_pattern(crystSys))
    #     temp = np.zeros((6, 6))
    #     temp_err = np.zeros((6, 6))
    #     for i in range(6):
    #         for j in range(6):
    #             Stress_string[i, j] = "C"+str(i+1)+str(j+1)
    #             temp[i, j], temp_err[i, j] = fit(i, j)
    #     temp = (temp+temp.T)/2
    #     temp_err = (temp_err+temp_err.T)/2
    #     for i in range(6):
    #         for j in range(i, 6):
    #             Cij_list[c[i, j]-1] = temp[i, j]
    #             Cij_errors_list[c[i, j]-1] = temp_err[i, j]
    #             index += 1

    # elif crystSys == 6:  # hexagonal
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[1, 0] = "C13"
    #     Stress_string[1, 1] = "C13"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 3] = "C44"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[6], Cij_errors_list[6] = fit(0, 1)  # c12
    #     fit02, error02 = fit(0, 2)  # 13
    #     fit10, error10 = fit(1, 0)  # 13
    #     fit11, error11 = fit(1, 1)  # 13
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     Cij_list[3], Cij_errors_list[3] = fit(1, 3)  # c44
    #     Cij_list[7] = (fit02+fit10+fit11)/3  # c13
    #     Cij_errors_list[7] = (error02+error10+error11)/3
    #     Cij_list[5] = 0.5*(Cij_list[0]-Cij_list[6])  # c66
    #     Cij_errors_list[5] = 0.5*(Cij_errors_list[0]-Cij_errors_list[6])

    # elif crystSys == 7:  # trignoal high
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C14"
    #     Stress_string[1, 0] = "C13 + C14"
    #     Stress_string[1, 1] = "C13 - C14"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 3] = "C44"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[6], Cij_errors_list[6] = fit(0, 1)  # c12
    #     fit02, error02 = fit(0, 2)  # 13
    #     fit03, error03 = fit(0, 3)  # 14
    #     fit10, error10 = fit(1, 0)  # 13+14
    #     fit11, error11 = fit(1, 1)  # 13-14
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     Cij_list[3], Cij_errors_list[3] = fit(1, 3)  # c44
    #     Cij_list[7] = (fit02+fit10+fit11)/3  # c13
    #     Cij_errors_list[7] = (error02+error10+error11)/3
    #     Cij_list[8] = (fit03+fit10-fit11)/3  # c14
    #     Cij_errors_list[8] = (error03+error10-error11)/3
    #     Cij_list[5] = 0.5*(Cij_list[0]-Cij_list[6])  # c66
    #     Cij_errors_list[5] = 0.5*(Cij_errors_list[0]-Cij_errors_list[6])

    # elif crystSys == 71:  # trignoal low
    #     Stress_string[0, 0] = "C11"
    #     Stress_string[0, 1] = "C12"
    #     Stress_string[0, 2] = "C13"
    #     Stress_string[0, 3] = "C14"
    #     Stress_string[0, 4] = "-C25"
    #     Stress_string[1, 0] = "C13 + C14"
    #     Stress_string[1, 1] = "C13 - C14"
    #     Stress_string[1, 2] = "C33"
    #     Stress_string[1, 3] = "C44"
    #     Stress_string[1, 5] = "C25"
    #     Cij_list[0], Cij_errors_list[0] = fit(0, 0)  # c11
    #     Cij_list[6], Cij_errors_list[6] = fit(0, 1)  # c12
    #     fit02, error02 = fit(0, 2)  # 13
    #     fit03, error03 = fit(0, 3)  # 14
    #     fit04, error04 = fit(0, 4)  # -25
    #     fit10, error10 = fit(1, 0)  # 13+14
    #     fit11, error11 = fit(1, 1)  # 13-14
    #     Cij_list[2], Cij_errors_list[2] = fit(1, 2)  # c33
    #     Cij_list[3], Cij_errors_list[3] = fit(1, 3)  # c44
    #     fit15, error15 = fit(1, 5)  # 25
    #     Cij_list[7] = (fit02+fit10+fit11)/3  # c13
    #     Cij_errors_list[7] = (error02+error10+error11)/3
    #     Cij_list[8] = (fit03+fit10-fit11)/3  # c14
    #     Cij_errors_list[8] = (error03+error10-error11)/3
    #     Cij_list[13] = (-fit04+fit15)/2  # c25
    #     Cij_errors_list[13] = (-error04+error15)/2
    #     Cij_list[5] = 0.5*(Cij_list[0]-Cij_list[6])  # c66
    #     Cij_errors_list[5] = 0.5*(Cij_errors_list[0]-Cij_errors_list[6])

    c, std_err = createCij()
    return c, std_err


def readCij():
    f = open("cij.dat", 'r')
    lines = f.readlines()
    f.close()
    c = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            c[i][j] = float(lines[i].strip().split()[j])
    return c


###################################################### Analysis ##########
class elast_consts:

    def __init__(self, arguments, cvoigt=np.zeros((6, 6)), std_err=np.zeros((6, 6))):
        print "\n----------------------------------------Analysis----------------------------------------"
        self.cvoigt = cvoigt
        self.std_err = std_err
        self.svoigt = np.linalg.inv(cvoigt)
        self.smat = self.getSmat()
        self.cmat = self.getCmat()
        if arguments.isPrintCijs:
            self.print_cvoigt()
        if arguments.isCheckStability:
            self.check_stability()
        if arguments.isCalcPolyModulus:
            self.calc_poly_modulus()
        if arguments.isWriteCijs:
            self.write_cvoigt()
        if arguments.isCalcMaxMin:
            self.minimum_elastic_moduli()
            self.maximum_elastic_moduli()
        if arguments.isCalcDirYoung:
            self.calc_dir_youngs_modulus()
        if arguments.isCalcDirLinCompress:
            self.calc_dir_lin_compress()
        self.show_dir_youngs_modulus(arguments.angles)

    def print_cvoigt(self):
        print "\nSymmetrized Elastic Constant (C +- err) (GPa):"
        for i in range(6):
            for j in range(6):
                print "%5.2f +- %5.2f\t" % (self.cvoigt[i, j], self.std_err[i, j]),
            print "\n",
        print ""

    def coeff(self, i, j):
        if i < 3 and j < 3:
            return 1.0
        if (i > 2 and j < 3) or (i < 3 and j > 2):
            return 0.5
        if i > 2 and j > 2:
            return 0.25

    def getSmat(self):
        smat = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        smat[i, j, k, l] = self.coeff(
                            voigt_mat[i, j], voigt_mat[k, l])*self.svoigt[voigt_mat[i, j], voigt_mat[k, l]]
        return smat

    def getCmat(self):
        voigt_mat = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
        cmat = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        cmat[i, j, k, l] = self.cvoigt[
                            voigt_mat[i, j], voigt_mat[k, l]]
        return cmat

    def check_stability(self):
        eigvals = np.linalg.eigvals(self.cvoigt)
        if any(eig < 0 for eig in eigvals):
            print "\nThe structure is instable! (Cij matrix has negative value)"
        else:
            print "\nThe structure is mechanically stable! (Cij matrix has no negative value)"
        print "Eigenvalues of Cij matrix:\n", eigvals

    def calc_poly_modulus(self):
        s = self.svoigt
        c = self.cvoigt
        Br = ((s[0, 0]+s[1, 1]+s[2, 2])+2*(s[0, 1]+s[0, 2]+s[1, 2]))**(-1)
        Bv = (1.0/9)*(c[0, 0]+c[1, 1]+c[2, 2]+2*(c[0, 1]+c[0, 2]+c[1, 2]))
        B = 0.5*(Bv+Br)
        Gr = 15*(4*(s[0, 0]+s[1, 1]+s[2, 2])-4 *
                 (s[0, 1]+s[0, 2]+s[1, 2])+3*(s[3, 3]+s[4, 4]+s[5, 5]))**(-1)
        Gv = (1.0/15)*((c[0, 0]+c[1, 1]+c[2, 2])+3 *
                       (c[3, 3]+c[4, 4]+c[5, 5])-(c[0, 1]+c[0, 2]+c[1, 2]))
        G = 0.5*(Gv+Gr)
        E = 9*B*G/(3*B+G)
        v = (3*B-2*G)/(6*B+2*G)
        print "\nIsotropic Elastic Modulus (GPa):"
        print "Young's Modulus: E=", "%6.2f" % E
        print "Shear Modulus: G=", "%6.2f" % G
        print "Bulk Modulus: B=", "%6.2f" % B
        print "Poisson's Ratio: v=", "%6.2f" % v

    # a and b are directional vectors
    def a(self, theta, phi):
        return np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)

    def b(self, theta, phi, chi):
        return np.cos(theta)*np.cos(chi)*np.cos(phi)-np.sin(chi)*np.sin(phi), np.cos(theta)*np.cos(chi)*np.sin(phi)+np.sin(chi)*np.cos(phi), -np.sin(theta)*np.cos(chi)

    def dir_youngs_moduli(self, angles):
        theta, phi = angles
        a = self.a(theta, phi)
        youngs_moduli = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        youngs_moduli += a[i]*a[j]*a[k] * \
                            a[l]*self.smat[i, j, k, l]
        return 1.0/youngs_moduli

    def show_dir_youngs_modulus(self, angles):
        if len(angles) != 0:
            angles = [float(ang) for ang in angles]
            angles_rad = [float(ang)*np.pi/180 for ang in angles]
            print "Young's Modulus along theta = %5.2f, phi = %5.2f : E = %6.2f" % (angles[0], angles[1], self.dir_youngs_moduli(angles_rad))

    def dir_shear_moduli(self, angles):
        theta, phi, chi = angles
        a = self.a(theta, phi)
        b = self.b(theta, phi, chi)
        shear_moduli = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        shear_moduli += a[i]*a[k]*b[j] * \
                            b[l]*self.smat[i, j, k, l]
        return 1.0/(4*shear_moduli)

    def minimum_elastic_moduli(self):
        from scipy.optimize import minimize
        x0 = np.array([0, 0, 0])
        x1 = np.array([0, 0])
        poisson_min = minimize(
            self.dir_poisson_ratio, x0, method='Powell', options={'xtol': 1e-8}).x
        shear_min = minimize(
            self.dir_shear_moduli, x0, method='Nelder-Mead', options={'xtol': 1e-8}).x
        youngs_min = minimize(
            self.dir_youngs_moduli, x1, method='Powell', options={'xtol': 1e-8}).x
        print "\nMinimum Modulus:\tEmin\tat\t(Theta,\tPhi)"
        print "Young's Modulus:\t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_youngs_moduli(youngs_min), youngs_min[0]*180/np.pi, youngs_min[1]*180/np.pi)
        print "Shear Modulus:  \t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_shear_moduli(shear_min), shear_min[0]*180/np.pi, shear_min[1]*180/np.pi)
        print "Poisson's Ratio:\t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_poisson_ratio(poisson_min), poisson_min[0]*180/np.pi, poisson_min[1]*180/np.pi)

    def maximum_elastic_moduli(self):
        from scipy.optimize import minimize
        dir_poisson_ratio = lambda angles: -self.dir_poisson_ratio(angles)
        dir_shear_moduli = lambda angles: -self.dir_shear_moduli(angles)
        dir_youngs_moduli = lambda angles: -self.dir_youngs_moduli(angles)
        x0 = np.array([0, 0, 0])
        x1 = np.array([0, 0])
        poisson_max = -1 * \
            minimize(
                dir_poisson_ratio, x0, method='Nelder-Mead', options={'xtol': 1e-8}).x
        shear_max = -1 * \
            minimize(
                dir_shear_moduli, x0, method='Nelder-Mead', options={'xtol': 1e-8}).x
        youngs_max = -1 * \
            minimize(
                dir_youngs_moduli, x1, method='Nelder-Mead', options={'xtol': 1e-8}).x
        print "\nMaximum Modulus:\tEmax\tat\t(Theta,\tPhi)"
        print "Young's Modulus:\t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_youngs_moduli(youngs_max), youngs_max[0]*180/np.pi, youngs_max[1]*180/np.pi)
        print "Shear Modulus:  \t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_shear_moduli(shear_max), shear_max[0]*180/np.pi, shear_max[1]*180/np.pi)
        print "Poisson's Ratio:\t%6.2f\tat\t%6.2f\t%6.2f" % (self.dir_poisson_ratio(poisson_max), poisson_max[0]*180/np.pi, poisson_max[1]*180/np.pi)

    def dir_poisson_ratio(self, angles):
        theta, phi, chi = angles
        a = self.a(theta, phi)
        b = self.b(theta, phi, chi)
        p_up = 0
        p_down = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        p_up += a[i]*a[j]*b[k]*b[l]*self.smat[i, j, k, l]
        return -p_up*self.dir_youngs_moduli([theta, phi])

    def dir_lin_compress(self, angles):
        theta, phi = angles
        a = self.a(theta, phi)
        lin_compress = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    lin_compress += a[i]*a[j]*self.smat[i, j, k, k]
        return lin_compress

    def calc_dir_youngs_modulus(self, npt=200):
        print "\nCalculating Directional Young's Modulus Data....."
        e = open("e.dat", 'w')
        p = np.linspace(0, 2*np.pi, npt)
        t = np.linspace(0, np.pi, npt)
        for theta in t:
            for phi in p:
                r = self.dir_youngs_moduli([theta, phi])
                e.write("%6.4f %6.4f %6.4f %6.4f\n" % (
                    r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta), r))
            e.write("\n")
        e.close()
        print "\nCalculating projections of Young's Modulus"
        p = np.linspace(0, 2*np.pi, npt)
        t = np.linspace(0, 2*np.pi, npt)
        np.savetxt('e_xy.dat', np.c_[
                   p, self.dir_youngs_moduli([np.pi/2, p])], delimiter='\t', fmt='%6.4f')
        np.savetxt('e_yz.dat', np.c_[
                   t, self.dir_youngs_moduli([t, np.pi/2])], delimiter='\t', fmt='%6.4f')
        np.savetxt('e_xz.dat', np.c_[
                   t, self.dir_youngs_moduli([t, 0])], delimiter='\t', fmt='%6.4f')
        print "Complete!"

    def calc_dir_lin_compress(self, npt=200):
        print "\nCalculating Directional Linear Compressibility Data....."
        beta = open("beta.dat", 'w')
        p = np.linspace(0, 2*np.pi, npt)
        t = np.linspace(0, np.pi, npt)
        for theta in t:
            for phi in p:
                r = self.dir_lin_compress([theta, phi])
                beta.write("%6.4f %6.4f %6.4f %6.4f\n" % (
                    r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta), r))
            beta.write("\n")
        beta.close()
        print "Complete!"

    def write_cvoigt(self):
        np.savetxt(
            'cij.dat', np.array(self.cvoigt), delimiter='\t', fmt='%5.2f')

    def write_svoigt(self):
        np.savetxt(
            'sij.dat', np.array(self.svoigt), delimiter='\t', fmt='%5.2f')
###################################################### Apply strain ######


def task_applystrain():
    deformation = get_deform_file()
    [p, bot, top] = readPoscar()
    shutil.move(poscar, 'POSCAR_backup')
    num_deform = len(deformation)
    for index_deform in range(num_deform):
        pos = apply_deform(p, deformation[index_deform])
        writePoscar(pos, bot, top, 'POSCAR'+str(index_deform+1))
        print "Applied deformation complete! (%d/%d)" % (index_deform+1, num_deform)
        index_deform = index_deform+1
    shutil.move('POSCAR_backup', poscar)


###################################################### Main Program ######


def main():
    options = argparse.ArgumentParser(
        description=info_text, formatter_class=argparse.RawTextHelpFormatter)
    options.add_argument(
        dest='task_opt', action='store', type=int, metavar='task_option',
        help="Select task (required):\n"
        "1)  Calculation preparation\n"
        "2)  Extract data from DFT calculation and analysis\n"
        "3)  Read Cijs from cij.dat file and anaylsis\n"
        "4)  Read strain from pattern.dat and apply strain")
    options.add_argument(
        dest='crystSys', action='store', type=int, metavar='crystal_system',
        help="Select crystal system (required):\n"
        "1)  Cubic\n"
        "2)  Tetragonal (4mm, -42m, 422, 4/mmm)\n"
        "21) Tetragonal (4, -4, 4/m)\n"
        "3)  Orthorhombic\n"
        "4)  Monoclinic (beta <> 90)\n"
        "41) Monoclinic (gamma <> 90)\n"
        "5)  Triclinic\n"
        "6)  Hexagonal\n"
        "7)  Trigonal (32, -3m, 3m)\n"
        "71) Trigonal (3, -3)")
    options.add_argument(
        '-n', dest='num_calcPoint', action='store', type=int, metavar='num_calcPoint', default=4,
        help="Number of calculation points per group of deformation (default: 4)")
    options.add_argument(
        '-e', dest='angles', action='append', metavar='theta_and_phi', default=[],
        help="Calculate Young's Modulus along specific direction (theta, phi) (first theta then phi in degree)")
    options.add_argument(
        '-d', dest='delta', action='store', type=float, metavar='delta', default=0.005,
        help="Magnitude of deformations intervals (default: 0.005 (0.5 percentage))")
    options.add_argument(
        '-cy', dest='isCalcDirYoung', action='store_true',
        help="Analysis: Calculate directional Young\'s modulus")
    options.add_argument(
        '-cl', dest='isCalcDirLinCompress', action='store_true',
        help="Analysis: Calculate directional linear compressiblity")
    options.add_argument(
        '-cm', dest='isCalcMaxMin', action='store_true',
        help="Analysis: Find the maximum and minimum of directional elastic modulus")
    options.add_argument(
        '-no-cs', dest='isCheckStability', action='store_false',
        help="Disable Analysis: Check stability")
    options.add_argument(
        '-no-cp', dest='isCalcPolyModulus', action='store_false',
        help="Disable Analysis: Calculate polycrystalline elastic modulus")
    options.add_argument(
        '-p', dest='isPrintCijs', action='store_false',
        help="Analysis: Print Cij")
    options.add_argument(
        '-no-cw', dest='isWriteCijs', action='store_false',
        help="Disable Analysis: Write Cij to cij.dat")
    options.add_argument(
        '-debug', dest='isDebug', action='store_false',
        help="Debug tag to write all strain and stress data to stress.txt")
    arguments = options.parse_args()
    printInfo()
    if arguments.task_opt == 1:
        try:
            elastic = preprocess(
                arguments.crystSys, arguments.num_calcPoint, arguments.delta)
        except:
            try:
                shutil.move('POSCAR_backup', poscar)
            except:
                pass
            sys.exit(1)
    elif arguments.task_opt == 2:
        postprocess(arguments)
    elif arguments.task_opt == 3:
        postprocess_read_cij(arguments)
    elif arguments.task_opt == 4:
        task_applystrain()
    elif arguments.task_opt == 5:
        postprocess_test(arguments)
    else:
        print "Error in Selection!"
        sys.exit(1)


# Invoke main program
if __name__ == '__main__':
    main()
