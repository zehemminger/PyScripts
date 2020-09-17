import csv
from array import array
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#PPM Error for fragments
PPM_Error_tolerance = 1

#N-Terminal Modification
N_Term_Mod = 0

#C-Terminal Mofidicaiton
C_Term_Mod = 0

#Atom Masses
n = 14.003074
o = 15.994915
h = 1.007825
p = 1.007276
car = 12.000000

#Import Amino Acid Sequence
Seq = open("Sequence.txt","r")
string = Seq.read()
List = list(string)
x = len(List)
length_sequence = len(List)
Seq.close()

#Import Ammino Acid masses
open("AminoAcids.csv","r")
AAs = []
MSs = []
with open("AminoAcids.csv") as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        AAs.append(row[0])
        MSs.append(row[1])
MSs_int= list(map(float, MSs))
#print(AAs)
#print(MSs)

print('calculating molecular weight')
#Calculating the molecular weight of the protein
Total_Mass = 0
t = 0
for i in range(0, x):
    while t != 20:
        if List[i] == AAs[t]:
            Total_Mass = Total_Mass + MSs_int[t]
            break
        else:
            t = t+1
    t = 0
Total_Mass_Two = Total_Mass + o + 2 * h
#print(Total_Mass_Two)
print("[M+H]+ = ", Total_Mass_Two + p)

print('generating termainal fragments')
#Generate a fragments
Fragment_Mass = 0
t = 0
arra = array('d')
arra_base = array('d')
arraf = []
aflength = 1
arrlenlistaf = []
stra = ''
for i in range(0, x):
    while t != 20:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != -1:
        stra = stra + List[i]
        i = i - 1
    arraf.append(stra)
    stra = ''
    Base_Fragment_Mass = Fragment_Mass - car - o
    Total_Fragment_Mass = Base_Fragment_Mass + p # Why need subtract H? Need to satisfy the carbonyl group. Maybe hooks to nitrogen?
    arra_base.append(Base_Fragment_Mass)
    arra.append(Total_Fragment_Mass)
    arrlenlistaf.append(aflength)
    aflength = aflength+1
    t = 0
#print("A_Fragments = ", arra)
#print("A_Fragments = ", arra_base)
#print("A_Numbers = ", arrlenlistaf)
#print("A_Fragments = ", arraf)
a = len(arra)

#Generate b fragments
Fragment_Mass = 0
t = 0
arrb = array('d')
arrb_base = array('d')
arrbf = []
bflength = 1
arrlenlistbf = []
strb = ''
for i in range(0, x):
    while t != 20:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != -1:
        strb = strb + List[i]
        i = i - 1
    arrbf.append(strb)
    strb = ''
    Base_Fragment_Mass = Fragment_Mass
    Total_Fragment_Mass = Base_Fragment_Mass + p # Why need subtract H? Need to satisfy the carbonyl group. Maybe hooks to nitrogen?
    arrb_base.append(Base_Fragment_Mass)
    arrb.append(Total_Fragment_Mass)
    arrlenlistbf.append(bflength)
    bflength = bflength+1
    t = 0
#print("B_Fragments = ", arrb)
#print("B_Fragments = ", arrb_base)
#print("B_Fragments = ", arrbf)
b = len(arrb)

#Generate c fragments
Fragment_Mass = 0
t = 0
arrc = array('d')
arrc_base = array('d')
arrcf = []
cflength = 1
arrlenlistcf = []
strc = ''
for i in range(0, x-1):
    while t != 20:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != -1:
        strc = strc + List[i]
        i = i - 1
    arrcf.append(strc)
    strc = ''
    Base_Fragment_Mass = Fragment_Mass + n + h + h + h
    Total_Fragment_Mass = Base_Fragment_Mass + p #Why add an extra H? Maybe when the bond breaks, another hydrogen is needed to satisfy the N that was separated. Same as y fragment
    arrc_base.append(Base_Fragment_Mass)
    arrc.append(Total_Fragment_Mass)
    arrlenlistcf.append(cflength)
    cflength = cflength+1
    t = 0
#print("C_Fragments = ", arrc)
#print("C_Fragments = ", arrc_base)
#print("C_Fragments = ", arrcf)
c = len(arrc)

#Generate x fragments
Fragment_Mass = 0
t = 0
arrx = array('d')
arrx_base = array('d')
arrxf = []
xflength = length_sequence
arrlenlistxf = []
strx = ''
for i in range(x-1, 0, -1):
    while t != 20:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != x:
        strx = strx + List[i]
        i = i + 1
    arrxf.append(strx)
    strx = ''
    Base_Fragment_Mass = Fragment_Mass + o + car + o
    Total_Fragment_Mass = Base_Fragment_Mass + p #Need to check this.
    arrx_base.append(Base_Fragment_Mass)
    arrx.append(Total_Fragment_Mass)
    arrlenlistxf.append(xflength)
    xflength = xflength - 1
    t = 0
#print("X_Fragments = ", arrx)
#print("X_Fragments = ", arrx_base)
#print("X_Fragments = ", arrxf)
x = len(arrx)

#Generate y fragments
Fragment_Mass = 0
t = 0
arry = array('d')
arry_base = array('d')
arryf = []
yflength = length_sequence
arrlenlistyf = []
stry = ''
for i in range(x, 0, -1):
    while t != 21:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != x+1:
        stry = stry + List[i]
        i = i + 1
    arryf.append(stry)
    stry = ''
    Base_Fragment_Mass = Fragment_Mass + o + h
    Total_Fragment_Mass = Base_Fragment_Mass + p + h #Why add an extra H? Maybe when the bond breaks, another hydrogen is needed to satisfy the N that was separated. Same as c fragment
    arry_base.append(Base_Fragment_Mass)
    arry.append(Total_Fragment_Mass)
    arrlenlistyf.append(yflength)
    yflength = yflength - 1
    t = 0
#print("Y_Fragments = ", arry)
#print("Y_Fragments = ", arry_base)
#print("Y_Fragments = ", arryf)
y = len(arry)

#Generate z fragments
Fragment_Mass = 0
t = 0
arrz = array('d')
arrz_base = array('d')
arrzf = []
zflength = length_sequence
arrlenlistzf = []
strz = ''
for i in range(x, -1, -1):
    while t != 20:
        if List[i] == AAs[t]:
            Fragment_Mass = Fragment_Mass + MSs_int[t]
            break
        else:
            t = t+1
    while i != x+1:
        strz = strz + List[i]
        i = i + 1
    arrzf.append(strz)
    strz = ''
    Base_Fragment_Mass = Fragment_Mass + o + h - n - h
    Total_Fragment_Mass = Base_Fragment_Mass + p #Correct
    arrz_base.append(Base_Fragment_Mass)
    arrz.append(Total_Fragment_Mass)
    arrlenlistzf.append(zflength)
    zflength = zflength - 1
    t = 0
#print("Z_Fragments = ", arrz)
#print("Z_Fragments = ", arrz_base)
#print("Z_Fragments = ", arrzf)
z = len(arrz)

print('generating internal fragments')
#Make internal fragments.
List = list(string)
w = 2
x = 0
v = 0
fl = 1
bh = len(List)
bl = len(List)
strintf = ''
arrintf = []
arrlenlistfl = []
arrlenlistbl = []
while v != b-1:
    fl = fl + 1
    del List[0]
    Listhold = List.copy()
    for i in range(x - b + w, x):
        bl = bl - 1
        List.reverse()
        del List[0]
        List.reverse()
        w = w + 1
        strintf = "".join(List)
        arrintf.append(strintf)
        arrlenlistfl.append(fl)
        arrlenlistbl.append(bl)
    bl = bh
    w = 2
    x = x + 1
    w = w + x
    v = v + 1
    List = Listhold.copy()
b = len(arrb)
List = list(string)
#print(arrlenlistfl)
#print(arrlenlistbl)
#print("intfraglength = ", intfraglength)

#a Fragments with x Fragments
axintfrag=0
nterm = 0
cterm = 0
arrax=array('d')
for i in range(0,a-1):
    for i in range(0,x-1):
        axintfrag = Total_Mass_Two - arra_base[nterm] - arrx_base[cterm] - 2 * h
        if nterm + cterm != b-2:
            axintfrag = axintfrag + p
            arrax.append(float(axintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
intfraglength = len(arrax)
#print("AX_Fragments = ", arrax)
#print("AX_Fragments = ", arrintf)

#a Fragments with y Fragments
ayintfrag=0
nterm = 0
cterm = 0
arrays=array('d')
for i in range(0,a-1):
    for i in range(0,y-1):
        ayintfrag = Total_Mass_Two - arra_base[nterm]-arry_base[cterm] - h
        if nterm + cterm != b-2:
            ayintfrag = ayintfrag + p
            arrays.append(float(ayintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("AY_Fragments = ", arrays)
#print("AY_Fragments = ", arrintf)

#a Fragments with z Fragments
azintfrag=0
nterm = 0
cterm = 0
arraz=array('d')
for i in range(0,c-1):
    for i in range(0,z-1):
        azintfrag = Total_Mass_Two - arra_base[nterm]-arrz_base[cterm] - 2*h
        if nterm + cterm != b-2:
            azintfrag = azintfrag + p
            arraz.append(azintfrag)
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("AZ_Fragments = ", arraz)
#print("AZ_Fragments = ", arrintf)

#b Fragments with x Fragments
bxintfrag=0
nterm = 0
cterm = 0
arrbx=array('d')
for i in range(0,b-1):
    for i in range(0,x-1):
        bxintfrag = Total_Mass_Two - arrb_base[nterm]- arrx_base[cterm] - 2*h
        if nterm + cterm != b-2:
            bxintfrag = bxintfrag + p
            arrbx.append(float(bxintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("BX_Fragments = ", arrbx)
#print("BX_Fragments = ", arrintf)

#b Fragments with y Fragments
byintfrag=0
nterm = 0
cterm = 0
arrby=array('d')
for i in range(0,b-1):
    for i in range(0,y-1):
        byintfrag = Total_Mass_Two - arrb_base[nterm] - arry_base[cterm] - h
        if nterm + cterm != b-2:
            byintfrag = byintfrag + p
            #print('{0:.6f}'.format(byintfrag))
            arrby.append(float(byintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("BY_Fragments = ", arrby)
#print("BY_Fragments = ", arrintf)

#b Fragments with z Fragments
bzintfrag=0
nterm = 0
cterm = 0
arrbz=array('d')
for i in range(0,b-1):
    for i in range(0,z-1):
        bzintfrag = Total_Mass_Two - arrb_base[nterm]-arrz_base[cterm]
        if nterm + cterm != b-2:
            bzintfrag = bzintfrag + p
            arrbz.append(float(bzintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("BZ_Fragments = ", arrbz)
#print("BZ_Fragments = ", arrintf)

#c Fragments with x Fragments
cxintfrag=0
nterm = 0
cterm = 0
arrcx=array('d')
for i in range(0,c-1):
    for i in range(0,x-1):
        cxintfrag = Total_Mass_Two - arrc_base[nterm]-arrx_base[cterm]
        if nterm + cterm != b-2:
            cxintfrag = cxintfrag + p
            arrcx.append(float(cxintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("CX_Fragments = ", arrcx)
#print("CX_Fragments = ", arrintf)

#c Fragments with y Fragments
cyintfrag=0
nterm = 0
cterm = 0
arrcy=array('d')
for i in range(0,c-1):
    for i in range(0,y-1):
        cyintfrag = Total_Mass_Two - arrc_base[nterm] - arry_base[cterm] + h
        if nterm + cterm != b-2:
            cyintfrag = cyintfrag + p
            arrcy.append(float(cyintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("CY_Fragments = ", arrcy)
#print("CY_Fragments = ", arrintf)

#c Fragments with z Fragments
czintfrag=0
nterm = 0
cterm = 0
arrcz=array('d')
for i in range(0,c-1):
    for i in range(0,z-1):
        czintfrag = Total_Mass_Two - arrc_base[nterm]-arrz_base[cterm] + h
        if nterm + cterm != b-2:
            czintfrag = czintfrag + p
            arrcz.append(float(czintfrag))
        else:
            break
        cterm = cterm + 1
    nterm = nterm + 1
    cterm = 0
#print("CZ_Fragments = ", arrcz)
#print("CZ_Fragments = ", arrintf)


# Import Instrument Measured Fragments
open("Fragments.csv","r")
instfrag = []
with open("Fragments.csv") as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        instfrag.append(row[0])
instfrags= list(map(float, instfrag))
#print("Measured Fragments =", instfrags)

# Import Modificaitons
open("Modifications.csv","r")
modlist = []
with open("Modifications.csv") as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        modlist.append(row[0])
modifications = list(map(float, modlist))
#print("Modificaitons = ", modifications)
ml = len(modifications)
#print(ml)

mod = 0
loop = 0

#Open a new file
with open('Output.csv', 'w') as csvFile:
    csvFile.close()

#Create a title row for matched fragemnts file
csvData = [['Frag_number', 'Frag Type', 'Modification', 'Term Mod', 'Observed Mass', 'Theoredical Mass', 'Start AA', 'End AA', 'Error', 'Sequence']]
with open('Output.csv', 'a') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(csvData)
csvFile.close()

#Create array of terminal fragment counts
list_of_counts_arr_terminal = array('f')
for i in range(0, length_sequence):
    list_of_counts_arr_terminal.append(0)
#print(list_of_counts_arr_terminal)

#Create array of internal fragment counts
list_of_counts_arr_internal = array('f')
for i in range(0, length_sequence):
    list_of_counts_arr_internal.append(0)
#print(list_of_counts_arr_internal)

frag_number = 0
print('start matching fragments')
#Start Matching and modification matching
for i in range(0, ml):
    mod = modifications[loop]
    loop = loop + 1
    #Match peaks with A fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, a):
            if abs((instfrags[g] - arra[m] - mod - N_Term_Mod)/(arra[m] + mod + N_Term_Mod)*1000000) < PPM_Error_tolerance:
                #print("A Fragment")
                #print(instfrags[g])
                #print(arra[m])
                Error = (instfrags[g] - arra[m] - mod - N_Term_Mod)/(arra[m] + mod + N_Term_Mod)*1000000
                #print(Error)
                #print(arraf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'A Fragment', str(mod),str(N_Term_Mod), str(instfrags[g]), str(arra[m]), str(1), str(arrlenlistaf[m]), str(Error), str(arraf[m])]]
                holdv = m
                for i in range(holdv, -1, -1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv - 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0

    #Match peaks with B fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, b):
            if abs((instfrags[g] - arrb[m] - mod - N_Term_Mod)/(arrb[m] + mod + N_Term_Mod)*1000000) < PPM_Error_tolerance:
                #print("B Fragment")
                #print(instfrags[g])
                #print(arrb[m])
                Error = (instfrags[g] - arrb[m] - mod - N_Term_Mod)/(arrb[m] + mod + N_Term_Mod)*1000000
                #print(Error)
                #print(arrbf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'B Fragment', str(mod), str(N_Term_Mod), str(instfrags[g]), str(arrb[m]), str(1), str(arrlenlistbf[m]), str(Error), str(arrbf[m])]]
                holdv = m
                for i in range(holdv, -1, -1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv - 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0

    #Match peaks with C fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, c):
            if abs((instfrags[g] - arrc[m] - mod - N_Term_Mod)/(arrc[m] + mod + N_Term_Mod)*1000000) < PPM_Error_tolerance:
                #print("C Fragment")
                #print(instfrags[g])
                #print(arrc[m])
                Error = (instfrags[g] - arrc[m] - mod - N_Term_Mod)/(arrc[m] + mod + N_Term_Mod)*1000000
                #print(Error)
                #print(arrcf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'C Fragment', str(mod), str(N_Term_Mod), str(instfrags[g]), str(arrc[m]), str(1), str(arrlenlistcf[m]), str(Error), str(arrcf[m])]]
                holdv = m
                for i in range(holdv, -1, -1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv - 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0


    #Match peaks with X fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, x):
            if abs((instfrags[g] - arrx[m] - mod - C_Term_Mod)/(arrx[m] + mod + C_Term_Mod)*1000000) < PPM_Error_tolerance:
                #print("X Fragment")
                #print(instfrags[g])
                #print(arrx[m])
                Error = (instfrags[g] - arrx[m] - mod - C_Term_Mod)/(arrx[m] + mod + C_Term_Mod)*1000000
                #print(Error)
                #print(arrxf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'X Fragment', str(mod), str(C_Term_Mod), str(instfrags[g]), str(arrx[m]), str(arrlenlistxf[m]), str(length_sequence), str(Error), str(arrxf[m])]]
                holdv = length_sequence - m - 1
                for i in range(holdv, length_sequence, 1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0

    #Match peaks with Y fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, y):
            if abs((instfrags[g] - arry[m] - mod - C_Term_Mod)/(arry[m] + mod + C_Term_Mod) * 1000000) < PPM_Error_tolerance:
                #print("Y Fragment")
                #print(instfrags[g])
                #print(arry[m])
                Error = (instfrags[g] - arry[m] - mod - C_Term_Mod)/(arry[m] + mod + C_Term_Mod) * 1000000
                #print(Error)
                #print(arryf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'Y Fragment', str(mod), str(C_Term_Mod), str(instfrags[g]), str(arry[m]), str(arrlenlistyf[m]), str(length_sequence), str(Error), str(arryf[m])]]
                holdv = length_sequence - m - 1
                for i in range(holdv, length_sequence, 1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0

    #Match peaks with Z fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, z):
            if abs((instfrags[g] - arrz[m] - mod - C_Term_Mod)/(arrz[m] + mod + C_Term_Mod) * 1000000) < PPM_Error_tolerance:
                #print("Z Fragment")
                #print(instfrags[g])
                #print(arrz[m])
                Error = (instfrags[g] - arrz[m] - mod - C_Term_Mod)/(arrz[m] + mod + C_Term_Mod) * 1000000
                #print(Error)
                #print(arrzf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'Z Fragment', str(mod), str(C_Term_Mod), str(instfrags[g]), str(arrz[m]), str(arrlenlistzf[m]), str(length_sequence), str(Error), str(arrzf[m])]]
                holdv = length_sequence - m - 1
                for i in range(holdv, length_sequence, 1):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0

    '''#Match peaks with AX fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrax[m] - mod)/(arrax[m] + mod)*1000000) < PPM_Error_tolerance:
                #print("AX Fragment")
                #print(instfrags[g])
                #print(arrax[m])
                Error = (instfrags[g] - arrax[m] - mod)/(arrax[m] + mod)*1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'AX Fragment', str(mod), str(0), str(instfrags[g]), str(arrax[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0'''

    '''#Match peaks with AY fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrays[m] - mod)/(arrays[m] + mod)*1000000) < PPM_Error_tolerance:
                #print("AY Fragment")
                #print(instfrags[g])
                #print(arrays[m])
                Error = (instfrags[g] - arrays[m] - mod)/(arrays[m] + mod)*1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'AY Fragment', str(mod), str(0), str(instfrags[g]), str(arrays[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0'''

    '''#Match peaks with AZ fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arraz[m] - mod)/(arraz[m] + mod)*1000000) < PPM_Error_tolerance:
                #print("AZ Fragment")
                #print(instfrags[g])
                #print(arraz[m])
                Error = (instfrags[g] - arraz[m] - mod)/(arraz[m] + mod)*1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'AZ Fragment', str(mod), str(0), str(instfrags[g]), str(arraz[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v=9
            m = m + 1
        g = g + 1
        m = 0'''

    '''# Match peaks with BX fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrbx[m] - mod) / (arrbx[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("BX Fragment")
                #print(instfrags[g])
                #print(arrbx[m])
                Error = (instfrags[g] - arrbx[m] - mod) / (arrbx[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'BX Fragment', str(mod), str(0), str(instfrags[g]), str(arrbx[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0'''

    # Match peaks with BY fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrby[m] - mod) / (arrby[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("BY Fragment")
                #print(instfrags[g])
                #print(arrby[m])
                Error = (instfrags[g] - arrby[m] - mod) / (arrby[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'BY Fragment', str(mod), str(0), str(instfrags[g]), str(arrby[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0

    '''# Match peaks with BZ fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrbz[m] - mod) / (arrbz[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("BZ Fragment")
                #print(instfrags[g])
                #print(arrbz[m])
                Error = (instfrags[g] - arrbz[m] - mod) / (arrbz[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'BZ Fragment', str(mod), str(0), str(instfrags[g]), str(arrbz[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0'''

    '''# Match peaks with CX fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrcx[m] - mod) / (arrcx[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("CX Fragment")
                #print(instfrags[g])
                #print(arrcx[m])
                Error = (instfrags[g] - arrcx[m] - mod) / (arrcx[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'CX Fragment', str(mod), str(0), str(instfrags[g]), str(arrcx[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0'''
    
    '''# Match peaks with CY fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrcy[m] - mod) / (arrcy[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("CY Fragment")
                #print(instfrags[g])
                #print(arrcy[m])
                Error = (instfrags[g] - arrcy[m] - mod) / (arrcy[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'CY Fragment', str(mod), str(0), str(instfrags[g]), str(arrcy[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), str(Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0'''

    '''# Match peaks with CZ fragments
    ifl = len(instfrags)
    g = 0
    m = 0
    for i in range(0, ifl):
        for i in range(0, intfraglength):
            if abs((instfrags[g] - arrcz[m] - mod) / (arrcz[m] + mod) * 1000000) < PPM_Error_tolerance:
                #print("CZ Fragment")
                #print(instfrags[g])
                #print(arrcz[m])
                Error = (instfrags[g] - arrcz[m] - mod) / (arrcz[m] + mod) * 1000000
                #print(Error)
                #print(arrintf[m])
                frag_number = frag_number + 1
                csvData = [[frag_number, 'CZ Fragment', str(mod), str(0), str(instfrags[g]), str(arrcz[m]), str(arrlenlistfl[m]), str(arrlenlistbl[m]), (Error), str(arrintf[m])]]
                holdv = arrlenlistfl[m] - 1
                holdvn = arrlenlistbl[m]
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] + 1
                    holdv = holdv + 1
                with open('Output.csv', 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerows(csvData)
                csvFile.close()
            else:
                v = 9
            m = m + 1
        g = g + 1
        m = 0'''
#print(list_of_counts_arr_terminal)
#print(list_of_counts_arr_internal)
print('done with matching')
#Delete Blank rows from Output.csv file and transfer to Matched_Fragments.csv file
df = pd.read_csv('Output.csv')
df.to_csv('Matched_Fragments.csv', index = False)
os.remove('Output.csv')

#Creat DataFrame from Matched_Fragments.csv
my_dataframe = pd.read_csv('Matched_Fragments.csv')
#print(my_dataframe)

#Count the number of fragments
count = len(my_dataframe)
#print(length_sequence)
'''#Delete Duplicate Rows and past duplicates in duplicates.csv
print('deleting duplicates')
row = 0
rowd = 0
for i in range(0, count):
    for i in range(0, count):
        if my_dataframe['Frag_number'][row] != my_dataframe['Frag_number'][rowd] and my_dataframe['Observed Mass'][row] == my_dataframe['Observed Mass'][rowd] and abs(my_dataframe['Error'][row]) < abs(my_dataframe['Error'][rowd]):
            if my_dataframe['Start AA'][rowd] == 1 or my_dataframe['End AA'][rowd] == length_sequence:
                holdv = int(my_dataframe['Start AA'][rowd]) - 1
                holdvn = int(my_dataframe['End AA'][rowd])
                for i in range(holdv, holdvn):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] - 1
                    holdv = holdv + 1
            else:
                holdv = int(my_dataframe['Start AA'][rowd]) - 1
                holdvn = int(my_dataframe['End AA'][rowd])
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] - 1
                    holdv = holdv + 1
            my_dataframe.iloc[my_dataframe.index[rowd]] = np.nan
            #print(my_dataframe)
        elif my_dataframe['Frag_number'][row] != my_dataframe['Frag_number'][rowd] and my_dataframe['Observed Mass'][row] == my_dataframe['Observed Mass'][rowd] and abs(my_dataframe['Error'][row]) > abs(my_dataframe['Error'][rowd]):
            if my_dataframe['Start AA'][row] == 1 or my_dataframe['End AA'][row] == length_sequence:
                holdv = int(my_dataframe['Start AA'][row]) - 1
                holdvn = int(my_dataframe['End AA'][row])
                for i in range(holdv, holdvn):
                    list_of_counts_arr_terminal[holdv] = list_of_counts_arr_terminal[holdv] - 1
                    holdv = holdv + 1
            else:
                holdv = int(my_dataframe['Start AA'][row]) - 1
                holdvn = int(my_dataframe['End AA'][row])
                for i in range(holdv, holdvn):
                    list_of_counts_arr_internal[holdv] = list_of_counts_arr_internal[holdv] - 1
                    holdv = holdv + 1
            my_dataframe.iloc[my_dataframe.index[row]] = np.nan
            #print(my_dataframe)
        else:
            nothing = 0
        rowd = rowd + 1
    rowd = 0
    row = row + 1
print(my_dataframe)
my_dataframe = my_dataframe[pd.notnull(my_dataframe['Start AA'])]
my_dataframe.to_csv('Matched_Fragments_V2.csv')'''

print('making figures')
#Add two count arrays
list_of_counts_arr_total = array('f')
for i in range(0, length_sequence):
    list_of_counts_arr_total.append(0)

#print(list_of_counts_arr_total)
for i in range (0, length_sequence):
    list_of_counts_arr_total[m] = list_of_counts_arr_terminal[m] + list_of_counts_arr_internal[m]
    m = m +1
#print(list_of_counts_arr_total)

#Open a new file to put counts in.
with open('Counts.csv', 'w') as csvFile:
    csvFile.close()
m = 0
counts = [['AA','Terminal', 'Internal','Total']]
with open('Counts.csv', 'a') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerows(counts)

for i in range(1, length_sequence+1):
    counts = [[List[m],str(list_of_counts_arr_terminal[m]),str(list_of_counts_arr_internal[m]),list_of_counts_arr_total[m]]]
    with open('Counts.csv', 'a') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(counts)
    m = m + 1

#Open a new file to put counts in.
with open('Counts.csv', 'w') as csvFile:
    csvFile.close()
m = 0
counts = [['AA','Terminal', 'Internal','Total']]
with open('Counts.csv', 'a') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerows(counts)

for i in range(1, length_sequence+1):
    counts = [[List[m],str(list_of_counts_arr_terminal[m]),str(list_of_counts_arr_internal[m]),list_of_counts_arr_total[m]]]
    with open('Counts.csv', 'a') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(counts)
    m = m + 1

df = pd.read_csv('Counts.csv')
df.to_csv('Count.csv', index = False)
os.remove('Counts.csv')

# Open data file
data = pd.read_csv('count.csv')

# Find max value
Total = data["Total"].tolist()

max = max(Total)

# Extract Amino Acid Sequence
AA_list = data["AA"].tolist()

i = len(AA_list)
while i % 25 != 0:
    if i % 25 == 0:
        break
    else:
        i = i +1
        AA_list.append(" ")

AA_table = pd.DataFrame(np.array(AA_list).reshape(-1,25),index=None,columns=None)

AA_table_np = AA_table.to_numpy()

# Extract Terminal Fragments
Term_Frag = data["Terminal"].tolist()

i = len(Term_Frag)
while i % 25 != 0:
    if i % 25 == 0:
        break
    else:
        i = i +1
        Term_Frag.append("0")

Term_FM = pd.to_numeric(Term_Frag)
Term_TM = pd.DataFrame(np.array(Term_FM).reshape(-1,25),index=None,columns=None)

Term_TM.to_csv('Term_TM.csv',index=None)

Term_corr = pd.read_csv('Term_TM.csv')

# Extract Internal Fragments
Int_Frag = data["Internal"].tolist()

i = len(Int_Frag)
while i % 25 != 0:
    if i % 25 == 0:
        break
    else:
        i = i +1
        Int_Frag.append("0")

Int_FM = pd.to_numeric(Int_Frag)
Int_TM = pd.DataFrame(np.array(Int_FM).reshape(-1,25),index=None,columns=None)

Int_TM.to_csv('Int_TM.csv',index=None)

Int_corr = pd.read_csv('Int_TM.csv')

# Plot heatmaps as figure
fig, (ax,ax2) = plt.subplots(nrows=2)
fig.subplots_adjust(hspace=0.2,)

sns.heatmap(Term_corr,
            cmap="Reds",
            ax=ax,
            vmin=0,
            vmax=None,
            square=True,
            linewidths=0.5,
            linecolor='white',
            cbar=False,
            xticklabels=False,
            yticklabels=False,
            annot=AA_table_np,
            fmt=''
            )

fig.colorbar(ax.collections[0],
             ax=ax,
             location="right",
             use_gridspec=False,
             pad=0.05,
             shrink=0.75,
             )

ax.set_title('Terminal Fragments',
             fontsize=10,
             loc='left'
             )

sns.heatmap(Int_corr,
            cmap="Blues",
            ax=ax2,
            vmin=0,
            vmax=None,
            square=True,
            linewidths=0.5,
            linecolor='white',
            cbar=False,
            xticklabels=False,
            yticklabels=False,
            annot=AA_table_np,
            fmt='')

fig.colorbar(ax2.collections[0],
             ax=ax2,
             location="right",
             use_gridspec=False,
             pad=0.05,
             shrink=0.75
             )

ax2.set_title('Internal Fragments',
              fontsize=10,
              loc='left'
              )

plt.show()