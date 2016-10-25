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

def cosgammag(lK, Rg, beta, delta, omegat):
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
    return (Rwx*lty-Rwy*ltx)/(lt*Rg)
    
def cosgammak(lK, Rg, beta, delta, omegat):
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
    return (sa*(Ray*ltz-Raz*lty) + ca*(Rax*lty-Ray*ltx))/lt
    
def vt(lK, Rg, beta, delta, omegat):
    alpha  = 0.5*np.pi-beta
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    ltx    =  ca*cotd + cb*lK - cot*Rg
    lty    =     sotd         - sot*Rg
    ltz    = -sa*cotd + sb*lK
    lt     =  np.sqrt(ltx*ltx+lty*lty+ltz*ltz)
    dltxdt = -ca*sotd         + sot*Rg
    dltydt =     cotd         - cot*Rg
    dltzdt =  sa*sotd
    return (ltx*dltxdt+lty*dltydt+ltz*dltzdt)/lt
    
f, (ax1, ax2) = plt.subplots(1, 2)

beta   = np.radians(30)
delta  = np.radians(0)
lK     = 2
Rg     = 1
omegat = np.arange(0, 360, 1)
ax1.plot(omegat,    cosgammak(lK, Rg, beta, delta, omegat), 'r') # gammak
ax1.plot(omegat, Rg*cosgammag(lK, Rg, beta, delta, omegat), 'b') # Rg*gammag
ax1.plot(omegat,    vt(lK, Rg, beta, delta, omegat), 'g')
ax1.plot(omegat,    cosgammak(lK, Rg, beta, delta, omegat) -
                 Rg*cosgammag(lK, Rg, beta, delta, omegat) -
                    vt(lK, Rg, beta, delta, omegat), 'k--')

ax1.set_xlim(0,360)
ax1.set_ylim(-1,1)
ax1.set_xticks(np.arange(0,361,90))
ax1.set_yticks(np.arange(-1,1.1,0.5))

delta  = np.radians(30)
ax2.plot(omegat,    cosgammak(lK, Rg, beta, delta, omegat), 'r') # gammak
ax2.plot(omegat, Rg*cosgammag(lK, Rg, beta, delta, omegat), 'b') # gammag
ax2.plot(omegat,    vt(lK, Rg, beta, delta, omegat), 'g')        # vt/(omega*R)

ax2.plot(omegat,    cosgammak(lK, Rg, beta, delta, omegat) -     # gammak
                 Rg*cosgammag(lK, Rg, beta, delta, omegat) -     # gammag
                    vt(lK, Rg, beta, delta, omegat), 'k--')      # vt/(omega*R)
    
ax2.set_xlim(0,360)
ax2.set_ylim(-1,1)
ax2.set_xticks(np.arange(0,361,90))
ax2.set_yticks(np.arange(-1,1.1,0.5))  
    
#
# Set axis labels
ax1.set_xlabel(r"Phase angle \$\omega t\, [^{\circ}]\$",labelpad=10)
ax1.set_ylabel(r"Dimensionless power \$[-]\$",labelpad=10)
ax2.set_xlabel(r"Phase angle \$\omega t\, [^{\circ}]\$",labelpad=10)
#
# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig14_diagram.svg")
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    