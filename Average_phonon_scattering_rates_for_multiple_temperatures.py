# Comprehensive folder path is the path that includes both main folders and secondary folders for different temperatures.

# Main folder must be named in the same way as the secondary folder (i.e, main folder name: T70K, secondary folder name: T70K).

# The characteristic length is by default = 1 nm; it can be changed from the same code.

# Temperature range is set up to 2000 K; it can be changed from the same code. 


from numpy import array, concatenate , atleast_1d
import numpy as np
import pandas as pd
from pathlib import Path
import os
from pylab import *
import pickle
import sys
import itertools
from colorama import init
from termcolor import cprint 
from pyfiglet import figlet_format
from PIL import Image, ImageFont, ImageDraw


###########################################################################################################################
# main folder
datadir =  input("Enter the comprehensive folder path: ")
results = []
###########################################################################################################################

for i in range(1,2000):
    # load the data
    files_1 = Path(datadir+'\T'+str(i)+'K')
    files_2 = Path(datadir+'\T'+str(i)+'K'+'\T'+str(i)+'K')

    if os.path.exists(files_1):
        a_0 = i
        a_0 = np.array([a_0])

        files_1= str(Path(files_1))
        files_2= str(Path(files_2))
        
        ####################################################################################################################
        # Required files:
        OMEGA = np.loadtxt(files_1+'\BTE.omega')   # angular frequency [rad/ps]
        Pi = np.loadtxt(files_2+'\BTE.w_3ph')  # phonon scattering rate [1/ps]
        Ni = np.loadtxt(files_2+'\BTE.w_3ph_normal')    # normal scattering rate [1/ps]
        Ui = np.loadtxt(files_2+'\BTE.w_3ph_Umklapp')   # umklapp scattering rate [1/ps]
        Ii = np.loadtxt(files_1+'\BTE.w_isotopic')    # isotope scattering rate [1/ps]
        Vi = np.loadtxt(files_1+'\BTE.v')    # group velocity [km/s]

        ####################################################################################################################
        rows = len(OMEGA)          # number of rows [OMEGA]
        columns = len(OMEGA[0])    # number of columns [OMEGA]

        # Group velocity:
        Vi_0 = Vi[:,0]
        Vi_0 = np.array(Vi_0)
        Vi_0.resize((columns, rows))
        df = pd.DataFrame(Vi_0)
        df_t = df.T
        Vi_0 = df_t

        Vi_1 = Vi[:,1]
        Vi_1 = np.array(Vi_1)
        Vi_1.resize((columns, rows))
        df = pd.DataFrame(Vi_1)
        df_t = df.T
        Vi_1 = df_t

        Vi_2 = Vi[:,2]
        Vi_2 = np.array(Vi_2)
        Vi_2.resize((columns, rows))
        df = pd.DataFrame(Vi_2)
        df_t = df.T
        Vi_2 = df_t

        Vi = np.concatenate((Vi_0, Vi_1, Vi_2), axis=1)

        Vi = Vi**2

        mag1 = (Vi[:,[0]]+Vi[:,[1]]+Vi[:,[2]])**0.5
        Vi[:,[0]] = mag1

        mag2 = (Vi[:,[3]]+Vi[:,[4]]+Vi[:,[5]])**0.5
        Vi[:,[1]] = mag2

        mag3 = (Vi[:,[6]]+Vi[:,[7]]+Vi[:,[8]])**0.5
        Vi[:,[2]] = mag3

        mag4 = (Vi[:,[9]]+Vi[:,[10]]+Vi[:,[11]])**0.5
        Vi[:,[3]] = mag4

        mag5 = (Vi[:,[12]]+Vi[:,[13]]+Vi[:,[14]])**0.5
        Vi[:,[4]] = mag5

        mag6 = (Vi[:,[15]]+Vi[:,[16]]+Vi[:,[17]])**0.5
        Vi[:,[5]] = mag6

        df = pd.DataFrame(Vi)
        df = df.drop(df.columns[[6,7,8,9,10,11,12,13,14,15,16,17]], axis = 1)
        Vi = df

        ####################################################################################################################
        OMEGA = OMEGA*1.e12    # convert from rad/ps to rad/s

        KB = 1.380649*1.e-23   # Boltzmann constant [J/K]  
        h_bar = 1.05457*1.e-34  # Planck's constant [J.s]
        T = int(i) # temperature [K]

        # Heat capacity
        e = np.exp(h_bar*OMEGA/(KB*T))   # exponential part 
        Ci = ((h_bar**2)*(OMEGA**2)/((T**2)*KB))*(e/(e-1)**2)    # heat capacity [J/K]

        df = pd.DataFrame(Ci)
        Ci = df.replace([np.NaN, -np.NaN], 0)  # replace NaN with 0


        # Sum of heat capacity
        Sum_Ci = Ci.sum().sum() # sum all values

        ####################################################################################################################
        # Phonon Scattering Rate (3ph)

        Pi = Pi[:,1]
        Pi = np.array(Pi)
        Pi.resize((columns, rows))
        df = pd.DataFrame(Pi)
        df_t = df.T
        Pi = df_t

        # Ci*Pi
        Ci_Pi = Ci*Pi

        # Sum of (Ci*Pi)
        Sum_Ci_Pi = Ci_Pi.sum().sum() # sum all values

        # Average Phonon Scattering Rate (APSR)
        APSR = Sum_Ci_Pi/Sum_Ci        

        a_1 = APSR
        a_1 = np.array([a_1])


        ####################################################################################################################
        # Lifetime (3ph) 

        Ti=1/Pi

        df = pd.DataFrame(Ti)
        Ti = df.replace([np.NaN, -np.NaN], 0)  # replace NaN with 0

        # Ci*Ti
        Ci_Ti = Ci*Ti

        # Sum of (Ci*Ti)
        Sum_Ci_Ti = Ci_Ti.sum().sum() # sum all values

        # Average Lifetime (ALT)
        ALT = Sum_Ci_Ti/Sum_Ci

        a_2 = ALT
        a_2 = np.array([a_2])


        ####################################################################################################################
        # Normal Scattering Rate (NSR)

        Ni = Ni[:,3]
        Ni = np.array(Ni)
        Ni.resize((columns, rows))
        df = pd.DataFrame(Ni)
        df_t = df.T
        Ni = df_t

        # Ci*Ni
        Ci_Ni = Ci*Ni

        # Sum of (Ci*Ni)
        Sum_Ci_Ni = Ci_Ni.sum().sum() # sum all values

        # Average Normal Scattering Rate (ANSR)
        ANSR = Sum_Ci_Ni/Sum_Ci

        a_3 = ANSR
        a_3 = np.array([a_3])


        ####################################################################################################################
        # Umklapp Scattering Rate (USR)

        Ui = Ui[:,3]
        Ui = np.array(Ui)
        Ui.resize((columns, rows))
        df = pd.DataFrame(Ui)
        df_t = df.T
        Ui = df_t

        # Ci*Ui
        Ci_Ui = Ci*Ui

        # Sum of (Ci*Ui)
        Sum_Ci_Ui = Ci_Ui.sum().sum() # sum all values

        # Average Umklapp Scattering Rate (AUSR)
        AUSR = Sum_Ci_Ui/Sum_Ci

        a_4 = AUSR
        a_4 = np.array([a_4])


        ####################################################################################################################
        # Isotope Scattering Rate (ISR)

        Ii = Ii[:,1]
        Ii = np.array(Ii)
        Ii.resize((columns, rows))
        df = pd.DataFrame(Ii)
        df_t = df.T
        Ii = df_t

        # Ci*Ii
        Ci_Ii = Ci*Ii

        # Sum of (Ci*Ii)
        Sum_Ci_Ii = Ci_Ii.sum().sum() # sum all values

        # Average Isotope Scattering Rate (AISR)
        AISR = Sum_Ci_Ii/Sum_Ci

        a_5 = AISR
        a_5 = np.array([a_5])


        ####################################################################################################################
        # Boundary Scattering Rate (ABSR) 

        Li = 1  # characteristic length [nm]
    

        Bi = (Vi/Li) # boundary scattering rate 

        # Ci*Bi
        Ci_Bi = Ci*Bi

        Sum_Ci_Bi = Ci_Bi.sum().sum() # sum all values


        # Average Boundary Scattering Rate (ABSR) 
        ABSR = Sum_Ci_Bi/Sum_Ci

        a_6 = ABSR
        a_6 = np.array([a_6])

        ####################################################################################################################
        file = []
        df = pd.DataFrame(file)
        file = [a_0, a_1, a_2, a_3, a_4, a_5, a_6]
        file = np.array(file)
        file = concatenate([atleast_1d(a) for a in [file]])
        df = pd.DataFrame(file)
        file = df
        rows_2 = len(file)          
        columns_2 = len(file[0]) 
        file = np.array(file)
        file.resize((1, rows_2))
        df = pd.DataFrame(file)
        file =np.array(file)
        results.append(file)

    else:
        i=i

results=np.concatenate((results))
df = pd.DataFrame(results)
results = df
df.columns =['T(K)', 'APSR_3p (1/ps)', 'ALT_3p (ps)', 'ANSR (1/ps)','AUSR(1/ps) ',  'AISR (1/ps)', 'ABSR (1/ps)']
df.to_excel(input("Enter the path where you would like to save the results in an excel file")+'/results.xlsx')

ShowText = 'DONE!'

font = ImageFont.truetype('arialbd.ttf', 15) 
size = font.getsize(ShowText)  
image = Image.new('1', size, 1) 
draw = ImageDraw.Draw(image)
draw.text((0, 0), ShowText, font=font) 
for rownum in range(size[1]): 
    line = []
    for colnum in range(size[0]):
        if image.getpixel((colnum, rownum)): line.append(' '),
        else: line.append('#'),
    print (''.join(line))
