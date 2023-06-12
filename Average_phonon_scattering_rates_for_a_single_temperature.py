# Average phonon scattering rates for a single temperature 
# Main folder path is the path that includes BTE.omega, BTE.w_isotopic, etc...
# Secondary folder path is the path that includes BTE.w_3ph, BTE.w_3ph_normal, etc... 

import numpy as np
import pandas as pd

# load the data:
datadir_1 = input("Enter the main folder path: ")
datadir_2 = input("Enter the secondary folder path: ")

###########################################################################################################################
# Required files:
OMEGA =np.loadtxt(datadir_1+'\BTE.omega')    # angular frequency [rad/ps]
Ii = np.loadtxt(datadir_1+'\BTE.w_isotopic')    # isotope scattering rate [1/ps]
Vi = np.loadtxt(datadir_1+'\BTE.v')    # group velocity [km/s]
Pi = np.loadtxt(datadir_2+'\BTE.w_3ph')    # phonon scattering rate [1/ps]
Ni = np.loadtxt(datadir_2+'\BTE.w_3ph_normal')    # normal scattering rate [1/ps]
Ui = np.loadtxt(datadir_2+'\BTE.w_3ph_Umklapp')   # umklapp scattering rate [1/ps]
###########################################################################################################################
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

############################################################################################################################
OMEGA = OMEGA*1.e12    # convert from rad/ps to rad/s

KB = 1.380649*1.e-23   # Boltzmann constant [J/K]  
h_bar = 1.05457*1.e-34  # Planck's constant [J.s]
T = float(input("Enter the temperature in [K]: ")) # temperature [K]

# Heat capacity
e = np.exp(h_bar*OMEGA/(KB*T))   # exponential part 
Ci = ((h_bar**2)*(OMEGA**2)/((T**2)*KB))*(e/(e-1)**2)    # heat capacity [J/K]

df = pd.DataFrame(Ci)
Ci = df.replace([np.NaN, -np.NaN], 0)  # replace NaN with 0


# Sum of heat capacity
Sum_Ci = Ci.sum().sum() # sum all values

###########################################################################################################################
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

print ('Average Phonon Scattering Rate (APSR) =', APSR)   # Print Average Phonon Scattering Rate 

###########################################################################################################################
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

print('Average Lifetime (ALT) =', ALT)   # Print Average Lifetime 

###########################################################################################################################
# Normal Scattering Rate (ANSR)

# Normal scattering:
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

print('Average Normal Scattering Rate (ANSR) =', ANSR)  # Print Average Normal Scattering Rate

###########################################################################################################################
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
print('Average Umklapp Scattering Rate (AUSR) =', AUSR)  # Print Average Umklapp Scattering Rate 

###########################################################################################################################
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

print('Average Isotope Scattering Rate (AISR) =', AISR)  # Print Average Isotope Scattering Rate (AISR)

###########################################################################################################################
# Boundary Scattering Rate (BSR) 

Li = float(input("Enter the characteristic length in [nm]: "))  # characteristic length [nm]
Bi = (Vi/Li) # boundary scattering rate for L

# Ci*Bi
Ci_Bi = Ci*Bi

Sum_Ci_Bi = Ci_Bi.sum().sum() # sum all values


# Average Boundary Scattering Rate (ABSR) 
ABSR = Sum_Ci_Bi/Sum_Ci

print('Average Boundary Scattering Rate (ABSR) at L =', ABSR)
