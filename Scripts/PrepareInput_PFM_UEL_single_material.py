#
# -*- coding: mbcs -*-
#
# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
#from abaqus import *
#from abaqusConstants import *
#--------------------------------------------------------------
import os
import sys
import fileinput
import math
import numpy as np


#Predefined function
def line_prepender_nf(filename,nfilename, line):
    with open(filename, 'r') as f:
        content = f.read()
        with open(nfilename,'a+') as f2:
            f2.write(line.rstrip('\n') + '\n' + content)

def is_not_blank(s):
    return bool(s and s.strip())

#from abaqus import getInput
#--------------------------------------------------------------
Beg = "*Heading\n" # Beginning of the input file
Para= "\n*Node" # Beginning of the input file
End = "*End Step" # End of the input file.

#Eletype = float(getInput('Eletype for fibre from 1 to 5 (C3D10=1;C3D8=2;C3D20R=3;CPE8R=4;CPE4=5)', '1'))
Eletype = 5.0

if (Eletype==1.0):
    txtStart = "*Element, type=C3D10" # Text before the element connectivity
elif (Eletype==2.0):
    txtStart = "*Element, type=C3D8" # Text before the element connectivity
elif (Eletype==3.0):
    txtStart = "*Element, type=C3D20R" # Text before the element connectivity
elif (Eletype==4.0):
    txtStart = "*Element, type=CPE8R" # Text before the element connectivity
elif (Eletype==5.0):
    txtStart = "*Element, type=CPE4" # Text before the element connectivity



txtEnd = "\n*Nset, nset=Part-1-1_Set-matrix, generate" # Text after the element connectivity

#INP = getInput('INP name:', 'TPB_Cell_2.inp')
INP = 'SEN_original_file.inp'# Original input file

#---------------------Fibre-----------------------------------------

UEL01 = '74000'
UEL02 = '0.30'
UEL03 = '0.19' #Xlc
UEL04 = '300' #Gc
UEL05 = '0'
UEL06 = '0'
UEL07 = '0'

#---------------------Matrix-----------------------------------------

UEL08 = '3500'
UEL09 = '0.35'
UEL10 = '0.001'
UEL11 = '0.0035' #Matrix Gc
UEL12 = '0'
UEL13 = '0'
UEL14 = '0'

#---------------------Homogenised Laminar-----------------------------------------

UEL15 = '11000'
UEL16 = '0.35'
UEL17 = '0.025'
UEL18 = '0.4' #Gc
UEL19 = '0'
UEL20 = '0'
UEL21 = '0'


# Remove file Job-2.inp (if it exists)
filenames = ['Job-2.txt','F0.txt', 'F1.txt', 'Elset_matrix.txt','Elset_fibre.txt','FibreElement.txt','MatrixElement.txt','LaminarElement.txt','Temp_FibreElement.txt','Temp_MatrixElement.txt','Temp_LaminarElement.txt','VisualElements','F2.txt','F3.txt']
for myfile in filenames:
    if os.path.isfile(myfile):
        os.remove(myfile)


# Generate the first file (F1)
f = open(INP)
content = f.read()
text = content.split(Beg)[1].split(Para)[0].strip()
final_file = open("F0.txt", "w")
final_file.write(text)
final_file.close()
f.close()

fileToSearch = 'F0.txt'
tempFile = open( fileToSearch, 'r+' )
for line in fileinput.input( fileToSearch ):
    tempFile.write(line.replace('** PART INSTANCE: Part-1-1', '*Parameter\nXlc=0.005*2\nGc=0.010\nSolution=0\n** PARTS'))
tempFile.close()

with open('F0.txt','a') as DataFile:
    DataFile.write(Para)

f = open(INP)
content = f.read()
text = content.split(Para)[1].split(txtStart)[0].strip()
final_file = open("F1.txt", "w")
final_file.write(text)
final_file.close()
f.close()

# Generate the second file (F2)
f = open(INP)
content = f.read()
text = content.split(txtStart)[1].split(txtEnd)[0].strip()
final_file = open("VisualElements.txt", "w")
final_file.write(text)
final_file.close()
f.close()

myDataFile = open('VisualElements.txt', 'r+') 
myDataLines = myDataFile.readlines() 
myPointsList = [eval(dataLine) for dataLine in myDataLines]

nelem = len(myPointsList)
myDataFile.close()

Order=math.floor(math.log10(nelem))
num0=10*10**Order+1

DataFile = open('F2.txt', 'a')
if (Eletype==1.0):  
    i = 0
    while i<nelem: 

            num=i+ nelem +1
            linea_partida = myDataLines[i]
            nodes = tuple(map(str,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            DataFile.writelines('%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s' % (num, a, b, c, d, e, f, g, h, j, k))

            i = i+1

elif (Eletype==2.0):  
    i = 0
    while i<nelem: 

            num=i+ nelem +1
            linea_partida = myDataLines[i]
            nodes = tuple(map(int,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if (i < 1):
                    DataFile.writelines('%d, %d, %d, %d, %d, %d, %d, %d, %d' % (num, a, b, c, d, e, f, g, h))
            elif (i > 0):
                    DataFile.writelines('\n%d, %d, %d, %d, %d, %d, %d, %d, %d' % (num, a, b, c, d, e, f, g, h))

            i = i+1

elif (Eletype==3.0):  
    i = 0
    while i<nelem: 

            num=i/2+nelem+1
            i1 = i+1
            linea_partida = myDataLines[i]
            nodes = tuple(map(str,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            j = nodes[9]
            k = nodes[10]
            a1 = nodes[11]
            b1 = nodes[12]
            c1 = nodes[13]
            d1 = nodes[14]
            e1 = nodes[15]
            linea_partida1 = myDataLines[i1]
            nodes1 = tuple(map(str,linea_partida1.split(",")))
            f1 = nodes1[0]
            g1 = nodes1[1]
            h1 = nodes1[2]
            j1 = nodes1[3]
            k1 = nodes1[4]
            DataFile.writelines('%d, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \n%s, %s, %s, %s, %s' % (num, a, b, c, d, e, f, g, h, j, k, a1, b1, c1, d1, e1, f1, g1, h1, j1, k1))

            i = i+2

elif (Eletype==4.0):  
    i = 0
    while i<nelem: 

            num=i+ nelem +1
            linea_partida = myDataLines[i]
            nodes = tuple(map(int,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            e = nodes[5]
            f = nodes[6]
            g = nodes[7]
            h = nodes[8]
            if (i < 1):
                    DataFile.writelines('%d, %d, %d, %d, %d, %d, %d, %d, %d' % (num, a, b, c, d, e, f, g, h))
            elif (i > 0):
                    DataFile.writelines('\n%d, %d, %d, %d, %d, %d, %d, %d, %d' % (num, a, b, c, d, e, f, g, h))

            i = i+1

elif (Eletype==5.0):  
    i = 0
    while i<nelem: 

            num=i+ nelem +1
            linea_partida = myDataLines[i]
            nodes = tuple(map(int,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            if (i < 1):
                    DataFile.writelines('%d, %d, %d, %d, %d' % (num, a, b, c, d))
            elif (i > 0):
                    DataFile.writelines('\n%d, %d, %d, %d, %d' % (num, a, b, c, d))

            i = i+1
DataFile.close()

myDataFile = open('VisualElements.txt', 'r+')
myDataLines = myDataFile.readlines()
myPointsList = [eval(dataLine) for dataLine in myDataLines]

nelem = len(myPointsList)

myDataFile.close()


# 4 Find the element set for the matrix
f = open(INP)
array=[]
txtStart3='*Elset, elset=Part-1-1_Set-matrix, generate'  # TO CHANGE
txtEnd3='** Section: Section-matrix' # TO CHANGE
content = f.read()
text = content.split(txtStart3)[1].split(txtEnd3)[0].strip()

array=(text.split(","))
array=[int(i) for i in array]
N=array[1]-array[0]+1
matrixN=np.linspace(array[0],array[1], num=N)
f.close()


DataFile = open('MatrixElement.txt', 'a')
i=0
while i<nelem:
    if (myPointsList[i][0] in matrixN): #Matching the element number
            linea_partida = myDataLines[i]
            nodes = tuple(map(int,linea_partida.split(",")))
            a = nodes[1]
            b = nodes[2]
            c = nodes[3]
            d = nodes[4]
            if (myPointsList[i][0] != matrixN[-1]):
                    DataFile.writelines('%d, %d, %d, %d, %d\n' % (i+1, a, b, c, d))
            elif (myPointsList[i][0] == matrixN[-1]):
                    DataFile.writelines('%d, %d, %d, %d, %d' % (i+1, a, b, c, d))
                    #print i
    i=i+1

DataFile.close()


#----------------------------------------------------------------------


# Generate the third file (F3)
f = open(INP)
content = f.read()
text = content.split(txtEnd)[1].split(End)[0].strip()
final_file = open("F3.txt", "w")
final_file.write(text)
final_file.close()
f.close()

# Change DEPVAR in F3
fileToSearch = 'F3.txt'
tempFile = open( fileToSearch, 'r+' )
for line in fileinput.input( fileToSearch ):
    tempFile.write(line.replace('*User Material, constants=1', '1, S11, S11\n2, S22, S22\n3, S33, S33\n4, S12, S12\n5, E11, E11\n6, E22, E22\n7, E33, E33\n8, E12, E12\n9, PHI, PHI\n10, H, H'))
tempFile.close()

#Replace the User-material definition
## When you define the user material properties, simply put 999. there
fileToSearch = 'F3.txt'
tempFile = open( fileToSearch, 'r+' )
for line in fileinput.input( fileToSearch ):
		tempFile.write(line.replace('999.,', '*User Material, constants=0'))
tempFile.close()

fileToSearch = 'F3.txt'
tempFile = open( fileToSearch, 'r+' )
for line in fileinput.input( fileToSearch ):
		tempFile.write(line.replace('*Output, field', '*Output, field, frequency=5'))
tempFile.close()



##Input the Matrix material properties
line='*User element, nodes=4, type=U1, properties=7, coordinates=2, var=40,\n  INTEGRATION=4, TENSOR=PSTRAIN\n1,2\n1,4\n*ELEMENT, TYPE=U1, ELSET=SOLID1\n'
line_prepender_nf('MatrixElement.txt','Temp_matrixElement.txt', line)
with open('Temp_matrixElement.txt','a') as DataFile:
    DataFile.writelines('\n*UEL PROPERTY, ELSET=SOLID1, MATERIAL=MATRIX\n')
    DataFile.writelines(UEL08)
    DataFile.writelines(', ')
    DataFile.writelines(UEL09)
    DataFile.writelines(', ')
    DataFile.writelines('<Xlc>')
    DataFile.writelines(', ')
    DataFile.writelines('<Gc>')
    DataFile.writelines(', ')
    DataFile.writelines(UEL12)
    DataFile.writelines(', ')
    DataFile.writelines(UEL13)
    DataFile.writelines(', ')
    DataFile.writelines('<Solution>')

##Input the virtual element properties
line='*Element, type=CPE4, elset=Visualization\n'
line_prepender_nf('F2.txt','Temp_F2.txt', line)
with open('Temp_F2.txt','a') as DataFile:
    DataFile.write(txtEnd)

filenames = ['F0.txt', 'F1.txt','Temp_MatrixElement.txt','Temp_F2.txt','F3.txt']
with open('Job-2.inp', 'w') as DataFile:
    DataFile.writelines(Beg)
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                DataFile.write(line)
        DataFile.writelines('\n')
    DataFile.writelines(End)

#
#Search should only happen once
fileToSearch = 'Job-2.inp'
tempFile = open(fileToSearch, 'r+')
for line in fileinput.input(fileToSearch):
    tempFile.write(line.replace('*Solid Section, elset=Part-1-1_Set-matrix, material=Material-1',
                                 '*Solid Section, elset=Visualization, material=Material-1'))
tempFile.close()


####Eliminate existing files###########################################
filenames2 = ['Temp_MatrixElement.txt','Temp_F2.txt']
for myfile in filenames2:
    if os.path.isfile(myfile):
        os.remove(myfile)

# print 'Combine_file_finished' #Python 2.7
print (Combine_file_finished) #Python 3.0