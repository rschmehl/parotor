import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import np
#
# Figure size
mpl.rcParams['figure.figsize'] = 11, 5.625
# Text and numbers are substituded by Inkscape's PDF+Latex export option.
# By varying the sizes we can adjust the positioning
# http://stackoverflow.com/questions/30201310/use-of-hyphen-or-minus-sign-in-matplotlib-versus-compatibility-with-latex
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams.update({'font.size': 15})
mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20)
#
# This is important for font substitution
mpl.rcParams['svg.fonttype'] = 'none'

def gammag(lK, Rg, beta, delta, omegat):
    alpha  = 0.5*np.pi-beta
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    ltx    =  np.cos(alpha)*cotd + cb*lK - cot*Rg
    lty    =                sotd         - sot*Rg
    ltz    = -np.sin(alpha)*cotd + sb*lK
    Rwx    =  cot*Rg
    Rwy    =  sot*Rg
    Rwz    =  0
    lt     =  np.sqrt(ltx*ltx+lty*lty+ltz*ltz)
    return np.degrees(np.arccos((Rwx*lty-Rwy*ltx)/(lt*Rg)))
    
def gammak(lK, Rg, beta, delta, omegat):
    alpha  = 0.5*np.pi-beta
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    Rax    =  ca*cotd
    Ray    =     sotd
    Raz    = -sa*cotd
    ltx    =  Rax + cb*lK - cot*Rg
    lty    =  Ray         - sot*Rg
    ltz    =  Raz + sb*lK
    lt     =  np.sqrt(ltx*ltx+lty*lty+ltz*ltz)
    return np.degrees(np.arccos((sa*(Ray*ltz-Raz*lty) + ca*(Rax*lty-Ray*ltx))/lt))
    
def gammak2(lK, Rg, beta, delta, omegat):
    alpha  = 0.5*np.pi-beta
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    cotdt  =  np.cos(np.radians(omegat+90)+delta)
    sotdt  =  np.sin(np.radians(omegat+90)+delta)    
    ltx    =  np.cos(alpha)*cotd + cb*lK - cot*Rg
    lty    =                sotd         - sot*Rg
    ltz    = -np.sin(alpha)*cotd + sb*lK
    eyrx   =  np.cos(alpha)*cotdt
    eyry   =                sotdt
    eyrz   = -np.sin(alpha)*cotdt
    lt     =  np.sqrt(ltx*ltx+lty*lty+ltz*ltz)
    return np.degrees(np.arccos((eyrx*ltx+eyry*lty+eyrz*ltz)/lt))
    
f, (ax1, ax2) = plt.subplots(1, 2)

lK     = 2
Rg     = 1
beta   = np.radians(30)
omegat = np.arange(0, 360, 1)

delta  = np.radians(0)
ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'r')
#ax1.plot(omegat, gammak(lK, Rg, beta, delta, omegat), 'b--')
#print "gammag @  90deg = ", gammag(lK, Rg, beta, delta,  90)
#print "gammag @ 270deg = ", gammag(lK, Rg, beta, delta, 270)
#print "gammag_max      = ", np.max(gammag(lK, Rg, beta, delta, omegat))
#print "gammag_min      = ", np.min(gammag(lK, Rg, beta, delta, omegat))
#print "max @          = ", omegat[np.argmax(gammag(lK, Rg, beta, delta, omegat))]
#print "min @          = ", omegat[np.argmin(gammag(lK, Rg, beta, delta, omegat))]
#print "gammag_max      = ", gammag(lK, Rg, beta, delta, omegat[np.argmax(gammag(lK, Rg, beta, delta, omegat))])
#print "gammag_min      = ", gammag(lK, Rg, beta, delta, omegat[np.argmin(gammag(lK, Rg, beta, delta, omegat))])
delta  = np.radians(30)
ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'r')
delta  = np.radians(60)
ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'r')
delta  = np.radians(90)
ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'r--')
delta  = np.radians(120)
#ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'b')
delta  = np.radians(150)
#ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'b')
delta  = np.radians(180)
#ax1.plot(omegat, gammag(lK, Rg, beta, delta, omegat), 'b--')

ax1.set_xlim(0,360)
ax1.set_ylim(0,180)
ax1.set_xticks(np.arange(0,361,90))
ax1.set_yticks(np.arange(0,181,30))

delta  = np.radians(0)
ax2.plot(omegat, gammak(lK, Rg, beta, delta, omegat), 'r')
#ax2.plot(omegat, gammak2(lK, Rg, beta, delta, omegat), 'b--')
#print "gammak @   0deg = ", gammak(lK, Rg, beta, delta, 0)
#print "gammak @  90deg = ", gammak(lK, Rg, beta, delta, 90)
#print "gammak @ 180deg = ", gammak(lK, Rg, beta, delta, 180)
#print "gammak @ 270deg = ", gammak(lK, Rg, beta, delta, 270)
delta  = np.radians(30)
ax2.plot(omegat, gammak(lK, Rg, beta, delta, omegat), 'r')
#ax2.plot(omegat, gammak2(lK, Rg, beta, delta, omegat), 'b--')
delta  = np.radians(60)
ax2.plot(omegat, gammak(lK, Rg, beta, delta, omegat), 'r')
delta  = np.radians(90)
ax2.plot(omegat, gammak(lK, Rg, beta, delta, omegat), 'r--')

ax2.set_xlim(0,360)
ax2.set_ylim(0,180)
ax2.set_xticks(np.arange(0,361,90))
ax2.set_yticks(np.arange(0,181,30))

#
# Set axis labels
ax1.set_xlabel(r"Phase angle \$\omega t\, [^{\circ}]\$",labelpad=10)
ax1.set_ylabel(r"Tether attachment angle \$\gammag\, [^{\circ}]\$",labelpad=10)
ax2.set_xlabel(r"Phase angle \$\omega t\, [^{\circ}]\$",labelpad=10)
#
# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig_torque_diagram.svg")
plt.show()
