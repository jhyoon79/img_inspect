import numpy as np

pc2m = 3.086e16

d_ngc4631 = 7 # Mpc
f_radio = 0.9816#1.2	# 0.9816	# Jy or 1.2 Jy from VLA
L_radio = f_radio*1e-26*4*np.pi*(d_ngc4631*1e6*pc2m)**2.

f_60um = 99.69# 85.4	Jy  or 99.69 Jy from IRAS
logL_60um = 6.014 + 2*np.log10(d_ngc4631)+np.log10(f_60um)
f_100um = 193.26	# 160.08	Jy or 193.26 Jy from IRAS

FIR = 1.26e-14*(2.58*f_60um+f_100um)
q = np.log10(FIR/3.75e12)-np.log10(f_radio*1e-26)

print np.log10(L_radio), logL_60um, q

erg2J = 1e-7
L_FIR = FIR/erg2J * 4*np.pi*(d_ngc4631*1e6*pc2m)**2.
SFR = 4.5e-44*L_FIR
log_SFRperRadio = np.log10(SFR/(f_radio*1e-26/erg2J*4*np.pi*(d_ngc4631*1e6*pc2m)**2.))+28

print SFR, log_SFRperRadio



