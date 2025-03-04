# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------------------
# get_chemical_symbol(sym/z) : IT RETURNS A SYMBOL, BASED ON THE DICTIONARY chemical_symbol[z]
# get_atomic_mass(sym)       : IT RETURNS THE ATOMIC MASS ON AMU, BASED IN atomic_mass[sym]
# get_covalent_radius(sym)   : IT RETURNS THE COVALENT RADIUS IN A, BASED IN covalent_radius[sym]
# check_atomic_sym(sym)      : IT RETURNS TRUE OF FALSE, BASED ON chemical_symbols
#------------------------------------------------------------------------------------------
#THE MAXIMUM ATOMIC NUMBER (Z) IS 86
#FOR CLARITY THE ATOMIC SYMBOL IT IS USED BY DEFAULT
chemical_symbol={}
chemical_symbol[0] ="X"
chemical_symbol[1] ="H"
chemical_symbol[2] ="He"
chemical_symbol[3] ="Li"
chemical_symbol[4] ="Be"
chemical_symbol[5] ="B"
chemical_symbol[6] ="C"
chemical_symbol[7] ="N"
chemical_symbol[8] ="O"
chemical_symbol[9] ="F"
chemical_symbol[10]="Ne"
chemical_symbol[11]="Na"
chemical_symbol[12]="Mg"
chemical_symbol[13]="Al"
chemical_symbol[14]="Si"
chemical_symbol[15]="P"
chemical_symbol[16]="S"
chemical_symbol[17]="Cl"
chemical_symbol[18]="Ar"
chemical_symbol[19]="K"
chemical_symbol[20]="Ca"
chemical_symbol[21]="Sc"
chemical_symbol[22]="Ti"
chemical_symbol[23]="V"
chemical_symbol[24]="Cr"
chemical_symbol[25]="Mn"
chemical_symbol[26]="Fe"
chemical_symbol[27]="Co"
chemical_symbol[28]="Ni"
chemical_symbol[29]="Cu"
chemical_symbol[30]="Zn"
chemical_symbol[31]="Ga"
chemical_symbol[32]="Ge"
chemical_symbol[33]="As"
chemical_symbol[34]="Se"
chemical_symbol[35]="Br"
chemical_symbol[36]="Kr"
chemical_symbol[37]="Rb"
chemical_symbol[38]="Sr"
chemical_symbol[39]="Y"
chemical_symbol[40]="Zr"
chemical_symbol[41]="Nb"
chemical_symbol[42]="Mo"
chemical_symbol[43]="Tc"
chemical_symbol[44]="Ru"
chemical_symbol[45]="Rh"
chemical_symbol[46]="Pd"
chemical_symbol[47]="Ag"
chemical_symbol[48]="Cd"
chemical_symbol[49]="In"
chemical_symbol[50]="Sn"
chemical_symbol[51]="Sb"
chemical_symbol[52]="Te"
chemical_symbol[53]="I"
chemical_symbol[54]="Xe"
chemical_symbol[55]="Cs"
chemical_symbol[56]="Ba"
chemical_symbol[57]="La"
chemical_symbol[58]="Ce"
chemical_symbol[59]="Pr"
chemical_symbol[60]="Nd"
chemical_symbol[61]="Pm"
chemical_symbol[62]="Sm"
chemical_symbol[63]="Eu"
chemical_symbol[64]="Gd"
chemical_symbol[65]="Tb"
chemical_symbol[66]="Dy"
chemical_symbol[67]="Ho"
chemical_symbol[68]="Er"
chemical_symbol[69]="Tm"
chemical_symbol[70]="Yb"
chemical_symbol[71]="Lu"
chemical_symbol[72]="Hf"
chemical_symbol[73]="Ta"
chemical_symbol[74]="W"
chemical_symbol[75]="Re"
chemical_symbol[76]="Os"
chemical_symbol[77]="Ir"
chemical_symbol[78]="Pt"
chemical_symbol[79]="Au"
chemical_symbol[80]="Hg"
chemical_symbol[81]="Tl"
chemical_symbol[82]="Pb"
chemical_symbol[83]="Bi"
chemical_symbol[84]="Po"
chemical_symbol[85]="At"
chemical_symbol[86]="Rn"

#------------------------------------------------------------------------------------------
###Units : Angstrom
covalent_radius={}
covalent_radius["H"] =0.32
covalent_radius["He"]=0.93
covalent_radius["Li"]=1.23
covalent_radius["Be"]=0.90
covalent_radius["B"] =0.82
covalent_radius["C"] =0.77
covalent_radius["N"] =0.75
covalent_radius["O"] =0.73
covalent_radius["F"] =0.72
covalent_radius["Ne"]=0.71
covalent_radius["Na"]=1.54
covalent_radius["Mg"]=1.36
covalent_radius["Al"]=1.18
covalent_radius["Si"]=1.11
covalent_radius["P"] =1.06
covalent_radius["S"] =1.02
covalent_radius["Cl"]=0.99
covalent_radius["Ar"]=0.98
covalent_radius["K"] =2.03
covalent_radius["Ca"]=1.74
covalent_radius["Sc"]=1.44
covalent_radius["Ti"]=1.32
covalent_radius["V"] =1.22
covalent_radius["Cr"]=1.18
covalent_radius["Mn"]=1.17
covalent_radius["Fe"]=1.17
covalent_radius["Co"]=1.16
covalent_radius["Ni"]=1.15
covalent_radius["Cu"]=1.17
covalent_radius["Zn"]=1.25
covalent_radius["Ga"]=1.26
covalent_radius["Ge"]=1.22
covalent_radius["As"]=1.20
covalent_radius["Se"]=1.16
covalent_radius["Br"]=1.14
covalent_radius["Kr"]=1.12
covalent_radius["Rb"]=2.16
covalent_radius["Sr"]=1.91
covalent_radius["Y"] =1.62
covalent_radius["Zr"]=1.45
covalent_radius["Nb"]=1.34
covalent_radius["Mo"]=1.30
covalent_radius["Tc"]=1.27
covalent_radius["Ru"]=1.25
covalent_radius["Rh"]=1.25
covalent_radius["Pd"]=1.28
covalent_radius["Ag"]=1.34
covalent_radius["Cd"]=1.48
covalent_radius["In"]=1.44
covalent_radius["Sn"]=1.41
covalent_radius["Sb"]=1.40
covalent_radius["Te"]=1.36
covalent_radius["I"] =1.33
covalent_radius["Xe"]=1.31
covalent_radius["Cs"]=2.35
covalent_radius["Ba"]=1.98
covalent_radius["La"]=1.69
covalent_radius["Ce"]=1.65
covalent_radius["Pr"]=1.65
covalent_radius["Nd"]=1.64
covalent_radius["Pm"]=1.63
covalent_radius["Sm"]=1.62
covalent_radius["Eu"]=1.85
covalent_radius["Gd"]=1.61
covalent_radius["Tb"]=1.59
covalent_radius["Dy"]=1.59
covalent_radius["Ho"]=1.58
covalent_radius["Er"]=1.57
covalent_radius["Tm"]=1.56
covalent_radius["Yb"]=1.74
covalent_radius["Lu"]=1.56
covalent_radius["Hf"]=1.44
covalent_radius["Ta"]=1.34
covalent_radius["W"] =1.30
covalent_radius["Re"]=1.28
covalent_radius["Os"]=1.26
covalent_radius["Ir"]=1.27
covalent_radius["Pt"]=1.30
covalent_radius["Au"]=1.34
covalent_radius["Hg"]=1.49
covalent_radius["Tl"]=1.48
covalent_radius["Pb"]=1.47
covalent_radius["Bi"]=1.46
covalent_radius["Po"]=1.46
covalent_radius["At"]=1.45
covalent_radius["Rn"]=1.90
#------------------------------------------------------------------------------------------
#Units: g/mol = amu
atomic_mass={}
atomic_mass["H"] = 1.0079
atomic_mass["He"]= 4.0026
atomic_mass["Li"]= 6.941
atomic_mass["Be"]= 9.0122
atomic_mass["B"] = 10.811
atomic_mass["C"] = 12.0107
atomic_mass["N"] = 14.0067
atomic_mass["O"] = 15.9994
atomic_mass["F"] = 18.9984
atomic_mass["Ne"]= 20.1797
atomic_mass["Na"]= 22.9897
atomic_mass["Mg"]= 24.305
atomic_mass["Al"]= 26.9815
atomic_mass["Si"]= 28.0855
atomic_mass["P"] = 30.9738
atomic_mass["S"] = 32.065
atomic_mass["Cl"]= 35.453
atomic_mass["Ar"]= 39.0983
atomic_mass["K"] = 39.948
atomic_mass["Ca"]= 40.078
atomic_mass["Sc"]= 44.9559
atomic_mass["Ti"]= 47.867
atomic_mass["V"] = 50.9415
atomic_mass["Cr"]= 51.9961
atomic_mass["Mn"]= 54.938
atomic_mass["Fe"]= 55.845
atomic_mass["Co"]= 58.6934
atomic_mass["Ni"]= 58.9332
atomic_mass["Cu"]= 63.546
atomic_mass["Zn"]= 65.39
atomic_mass["Ga"]= 69.723
atomic_mass["Ge"]= 72.64
atomic_mass["As"]= 74.9216
atomic_mass["Se"]= 78.96
atomic_mass["Br"]= 79.904
atomic_mass["Kr"]= 83.8
atomic_mass["Rb"]= 85.4678
atomic_mass["Sr"]= 87.62
atomic_mass["Y"] = 88.9059
atomic_mass["Zr"]= 91.224
atomic_mass["Nb"]= 92.9064
atomic_mass["Mo"]= 95.94
atomic_mass["Tc"]= 98.0
atomic_mass["Ru"]= 101.07
atomic_mass["Rh"]= 102.9055
atomic_mass["Pd"]= 106.42
atomic_mass["Ag"]= 107.8682
atomic_mass["Cd"]= 112.411
atomic_mass["In"]= 114.818
atomic_mass["Sn"]= 118.71
atomic_mass["Sb"]= 121.76
atomic_mass["Te"]= 126.9045
atomic_mass["I"] = 127.6
atomic_mass["Xe"]= 131.293
atomic_mass["Cs"]= 132.9055
atomic_mass["Ba"]= 137.327
atomic_mass["La"]= 138.9055
atomic_mass["Ce"]= 140.116
atomic_mass["Pr"]= 140.9077
atomic_mass["Nd"]= 144.24
atomic_mass["Pm"]= 145.0
atomic_mass["Sm"]= 150.36
atomic_mass["Eu"]= 151.964
atomic_mass["Gd"]= 157.25
atomic_mass["Tb"]= 158.9253
atomic_mass["Dy"]= 162.5
atomic_mass["Ho"]= 164.9303
atomic_mass["Er"]= 167.259
atomic_mass["Tm"]= 168.9342
atomic_mass["Yb"]= 173.04
atomic_mass["Lu"]= 174.967
atomic_mass["Hf"]= 178.49
atomic_mass["Ta"]= 180.9479
atomic_mass["W"] = 183.84
atomic_mass["Re"]= 186.207
atomic_mass["Os"]= 190.23
atomic_mass["Ir"]= 192.217
atomic_mass["Pt"]= 195.078
atomic_mass["Au"]= 196.9665
atomic_mass["Hg"]= 200.59
atomic_mass["Tl"]= 204.3833
atomic_mass["Pb"]= 207.2
atomic_mass["Bi"]= 208.9804
atomic_mass["Po"]= 209.0
atomic_mass["At"]= 210.0
atomic_mass["Rn"]= 222.0
#------------------------------------------------------------------------------------------    

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
#------------------------------------------------------------------------------------------
""" ROBUSTS, sym COULD BE a SYMBOL ("H") OR an INTEGER (1)
    HOW TO USE:
    >>> from utils.atomic import get_chemical_symbol
    """

def get_chemical_symbol(atomic_number):

    symbol = [ ]
    for i in atomic_number:
        if is_number(i):
            isymbol = chemical_symbol[int(i)]
            symbol.append(isymbol)
        else:
            isymbol = i
            symbol.append(isymbol)
        
    return symbol

#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from utils.atomic import get_atomic_mass
    """

def get_atomic_mass(sym):
    atomicmass=atomic_mass[sym]
    return atomicmass
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from utils.atomic import get_covalent_radius
    """

def get_covalent_radius(sym):
    covalentradius=covalent_radius[sym]
    return covalentradius
#------------------------------------------------------------------------------------------
chemical_symbols = [
'X',
'H',  #1  
'He', #2
'Li', #3
'Be', #4
'B',  #5
'C',  #6
'N',  #7
'O',  #8
'F',  #9
'Ne', #10
'Na', #11
'Mg', #12
'Al', #13
'Si', #14
'P',  #15
'S',  #16
'Cl', #17
'Ar', #18
'K',  #19
'Ca', #20
'Sc', #21
'Ti', #22
'V',  #23
'Cr', #24
'Mn', #25
'Fe', #26
'Co', #27
'Ni', #28
'Cu', #29
'Zn', #30
'Ga', #31
'Ge', #32
'As', #33
'Se', #34
'Br', #35 
'Kr', #36
'Rb', #37
'Sr', #38
'Y',  #39
'Zr', #40
'Nb', #41
'Mo', #42
'Tc', #43
'Ru', #44
'Rh', #45
'Pd', #46
'Ag', #47
'Cd', #48
'In', #49
'Sn', #50
'Sb', #51
'Te', #52
'I',  #53
'Xe', #54
'Cs', #55
'Ba', #56
'La', #57
'Ce', #58
'Pr', #59
'Nd', #60
'Pm', #61
'Sm', #62
'Eu', #63
'Gd', #64
'Tb', #65
'Dy', #66
'Ho', #67
'Er', #68
'Tm', #69
'Yb', #70
'Lu', #71
'Hf', #72
'Ta', #73
'W',  #74
'Re', #75
'Os', #76
'Ir', #77
'Pt', #78
'Au', #79
'Hg', #80
'Tl', #81
'Pb', #82
'Bi', #83
'Po', #84
'At', #85
'Rn',] #86
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from utils.atomic import check_atomic_sym
    """

def check_atomic_sym(sym):
    val0=False
    for ii in range(1,len(chemical_symbols)):
        if sym == chemical_symbols[ii]: val0=True
    return val0
#------------------------------------------------------------------------------------------
sym2atomic={}
sym2atomic['H'] =1
sym2atomic['He']=2 
sym2atomic['Li']=3 
sym2atomic['Be']=4 
sym2atomic['B'] =5
sym2atomic['C'] =6
sym2atomic['N'] =7
sym2atomic['O'] =8
sym2atomic['F'] =9
sym2atomic['Ne']=10 
sym2atomic['Na']=11 
sym2atomic['Mg']=12 
sym2atomic['Al']=13 
sym2atomic['Si']=14 
sym2atomic['P'] =15
sym2atomic['S'] =16
sym2atomic['Cl']=17 
sym2atomic['Ar']=18 
sym2atomic['K'] =19
sym2atomic['Ca']=20 
sym2atomic['Sc']=21 
sym2atomic['Ti']=22 
sym2atomic['V'] =23
sym2atomic['Cr']=24 
sym2atomic['Mn']=25 
sym2atomic['Fe']=26 
sym2atomic['Co']=27 
sym2atomic['Ni']=28 
sym2atomic['Cu']=29 
sym2atomic['Zn']=30 
sym2atomic['Ga']=31 
sym2atomic['Ge']=32 
sym2atomic['As']=33 
sym2atomic['Se']=34 
sym2atomic['Br']=35 
sym2atomic['Kr']=36 
sym2atomic['Rb']=37 
sym2atomic['Sr']=38 
sym2atomic['Y'] =39
sym2atomic['Zr']=40 
sym2atomic['Nb']=41 
sym2atomic['Mo']=42 
sym2atomic['Tc']=43 
sym2atomic['Ru']=44 
sym2atomic['Rh']=45 
sym2atomic['Pd']=46 
sym2atomic['Ag']=47 
sym2atomic['Cd']=48 
sym2atomic['In']=49 
sym2atomic['Sn']=50 
sym2atomic['Sb']=51 
sym2atomic['Te']=52 
sym2atomic['I'] =53
sym2atomic['Xe']=54 
sym2atomic['Cs']=55 
sym2atomic['Ba']=56 
sym2atomic['La']=57 
sym2atomic['Ce']=58 
sym2atomic['Pr']=59 
sym2atomic['Nd']=60 
sym2atomic['Pm']=61 
sym2atomic['Sm']=62 
sym2atomic['Eu']=63 
sym2atomic['Gd']=64 
sym2atomic['Tb']=65 
sym2atomic['Dy']=66 
sym2atomic['Ho']=67 
sym2atomic['Er']=68 
sym2atomic['Tm']=69 
sym2atomic['Yb']=70 
sym2atomic['Lu']=71 
sym2atomic['Hf']=72 
sym2atomic['Ta']=73 
sym2atomic['W'] =74
sym2atomic['Re']=75 
sym2atomic['Os']=76 
sym2atomic['Ir']=77 
sym2atomic['Pt']=78 
sym2atomic['Au']=79 
sym2atomic['Hg']=80 
sym2atomic['Tl']=81 
sym2atomic['Pb']=82 
sym2atomic['Bi']=83 
sym2atomic['Po']=84 
sym2atomic['At']=85 
sym2atomic['Rn']=86
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from utils.atomic import get_atomic_number
    """
 
def get_atomic_number(sym):
    atomicnumber=sym2atomic[sym]
    return atomicnumber
#------------------------------------------------------------------------------------------
#THE MAXIMUM ATOMIC NUMBER (Z) IS 86
#FOR CLARITY THE ATOMIC SYMBOL IT IS USED BY DEFAULT
chemical_color={}
chemical_color['H'] ='<1.00000000,1.00000000,1.00000000>'
chemical_color['He']='<0.85098039,1.00000000,1.00000000>'
chemical_color['Li']='<0.80000000,0.50196078,1.00000000>'
chemical_color['Be']='<0.76078431,1.00000000,0.00000000>'
chemical_color['B'] ='<1.00000000,0.70980392,0.70980392>'
chemical_color['C'] ='<0.56470588,0.56470588,0.56470588>'
chemical_color['N'] ='<0.18823529,0.31372549,0.97254902>'
chemical_color['O'] ='<1.00000000,0.05098039,0.05098039>'
chemical_color['F'] ='<0.56470588,0.87843137,0.31372549>'
chemical_color['Ne']='<0.70196078,0.89019608,0.96078431>'
chemical_color['Na']='<0.67058824,0.36078431,0.94901961>'
chemical_color['Mg']='<0.54117647,1.00000000,0.00000000>'
chemical_color['Al']='<0.74901961,0.65098039,0.65098039>'
chemical_color['Si']='<0.94117647,0.78431373,0.62745098>'
chemical_color['P'] ='<1.00000000,0.50196078,0.00000000>'
chemical_color['S'] ='<1.00000000,1.00000000,0.18823529>'
chemical_color['Cl']='<0.12156863,0.94117647,0.12156863>'
chemical_color['Ar']='<0.50196078,0.81960784,0.89019608>'
chemical_color['K'] ='<0.56078431,0.25098039,0.83137255>'
chemical_color['Ca']='<0.23921569,1.00000000,0.00000000>'
chemical_color['Sc']='<0.90196078,0.90196078,0.90196078>'
##chemical_color['Ti']='<0.74901961,0.76078431,0.78039216>'
chemical_color['Ti']='<0.45100000,0.76100000,0.98400000>'
chemical_color['V'] ='<0.65098039,0.65098039,0.67058824>'
chemical_color['Cr']='<0.54117647,0.60000000,0.78039216>'
chemical_color['Mn']='<0.61176471,0.47843137,0.78039216>'
chemical_color['Fe']='<0.87843137,0.40000000,0.20000000>'
chemical_color['Co']='<0.94117647,0.56470588,0.62745098>'
chemical_color['Ni']='<0.31372549,0.81568627,0.31372549>'
chemical_color['Cu']='<0.78431373,0.50196078,0.20000000>'
chemical_color['Zn']='<0.49019608,0.50196078,0.69019608>'
chemical_color['Ga']='<0.76078431,0.56078431,0.56078431>'
chemical_color['Ge']='<0.40000000,0.56078431,0.56078431>'
chemical_color['As']='<0.74117647,0.50196078,0.89019608>'
chemical_color['Se']='<1.00000000,0.63137255,0.00000000>'
chemical_color['Br']='<0.65098039,0.16078431,0.16078431>'
chemical_color['Kr']='<0.36078431,0.72156863,0.81960784>'
chemical_color['Rb']='<0.43921569,0.18039216,0.69019608>'
chemical_color['Sr']='<0.00000000,1.00000000,0.00000000>'
chemical_color['Y'] ='<0.58039216,1.00000000,1.00000000>'
chemical_color['Zr']='<0.58039216,0.87843137,0.87843137>'
chemical_color['Nb']='<0.45098039,0.76078431,0.78823529>'
chemical_color['Mo']='<0.32941176,0.70980392,0.70980392>'
chemical_color['Tc']='<0.23137255,0.61960784,0.61960784>'
chemical_color['Ru']='<0.14117647,0.56078431,0.56078431>'
chemical_color['Rh']='<0.03921569,0.49019608,0.54901961>'
chemical_color['Pd']='<0.00000000,0.41176471,0.52156863>'
chemical_color['Ag']='<0.75294118,0.75294118,0.75294118>'
chemical_color['Cd']='<1.00000000,0.85098039,0.56078431>'
chemical_color['In']='<0.65098039,0.45882353,0.45098039>'
chemical_color['Sn']='<0.40000000,0.50196078,0.50196078>'
chemical_color['Sb']='<0.61960784,0.38823529,0.70980392>'
chemical_color['Te']='<0.83137255,0.47843137,0.00000000>'
chemical_color['I'] ='<0.58039216,0.00000000,0.58039216>'
chemical_color['Xe']='<0.25882353,0.61960784,0.69019608>'
chemical_color['Cs']='<0.34117647,0.09019608,0.56078431>'
chemical_color['Ba']='<0.00000000,0.78823529,0.00000000>'
chemical_color['La']='<0.43921569,0.83137255,1.00000000>'
chemical_color['Ce']='<1.00000000,1.00000000,0.78039216>'
chemical_color['Pr']='<0.85098039,1.00000000,0.78039216>'
chemical_color['Nd']='<0.78039216,1.00000000,0.78039216>'
chemical_color['Pm']='<0.63921569,1.00000000,0.78039216>'
chemical_color['Sm']='<0.56078431,1.00000000,0.78039216>'
chemical_color['Eu']='<0.38039216,1.00000000,0.78039216>'
chemical_color['Gd']='<0.27058824,1.00000000,0.78039216>'
chemical_color['Tb']='<0.18823529,1.00000000,0.78039216>'
chemical_color['Dy']='<0.12156863,1.00000000,0.78039216>'
chemical_color['Ho']='<0.00000000,1.00000000,0.61176471>'
chemical_color['Er']='<0.00000000,0.90196078,0.45882353>'
chemical_color['Tm']='<0.00000000,0.83137255,0.32156863>'
chemical_color['Yb']='<0.00000000,0.74901961,0.21960784>'
chemical_color['Lu']='<0.00000000,0.67058824,0.14117647>'
chemical_color['Hf']='<0.30196078,0.76078431,1.00000000>'
chemical_color['Ta']='<0.30196078,0.65098039,1.00000000>'
chemical_color['W'] ='<0.12941176,0.58039216,0.83921569>'
chemical_color['Re']='<0.14901961,0.49019608,0.67058824>'
chemical_color['Os']='<0.14901961,0.40000000,0.58823529>'
chemical_color['Ir']='<0.09019608,0.32941176,0.52941176>'
chemical_color['Pt']='<0.81568627,0.81568627,0.87843137>'
chemical_color['Au']='<1.00000000,0.81960784,0.13725490>'
chemical_color['Hg']='<0.72156863,0.72156863,0.81568627>'
chemical_color['Tl']='<0.65098039,0.32941176,0.30196078>'
chemical_color['Pb']='<0.34117647,0.34901961,0.38039216>'
chemical_color['Bi']='<0.61960784,0.30980392,0.70980392>'
chemical_color['Po']='<0.67058824,0.36078431,0.00000000>'
chemical_color['At']='<0.45882353,0.30980392,0.27058824>'
chemical_color['Rn']='<0.25882353,0.50980392,0.58823529>'
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from utils.atomic import get_chemical_color
    """

def get_chemical_color(sym):
    color=chemical_color[sym]
    return color
#------------------------------------------------------------------------------------------
