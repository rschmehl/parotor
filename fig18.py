import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import np
from matplotlib import cm
from scipy.optimize import fmin
from matplotlib.ticker import MultipleLocator
import scipy.ndimage
import numpy.ma as ma
#
# Figure size
mpl.rcParams['figure.figsize'] = 10, 5.625
# Text and numbers are substituded by Inkscape's PDF+Latex export option.
# By varying the sizes we can adjust the positioning
# http://stackoverflow.com/questions/30201310/use-of-hyphen-or-minus-sign-in-matplotlib-versus-compatibility-with-latex
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams.update({'font.size': 22})
mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20)
#
# This is important for font substitution
mpl.rcParams['svg.fonttype'] = 'none'

class FormatFaker(object):
    def __init__(self, str): self.str = str
    def __mod__(self, stuff): return self.str
  
def Ftrx(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #
    ltx    =  ekxx*cotd + ekyx*sotd + cb*lK - cot*Rg
    lty    =  ekxy*cotd + ekyy*sotd         - sot*Rg
    ltz    =  ekxz*cotd             + sb*lK
    #
    return (ltx*ekxx+lty*ekxy+ltz*ekxz)/(ltx*ekzx+lty*ekzy+ltz*ekzz)

def Ftrxsum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Ftrx(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Ftrx(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Ftrx(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Ftrx(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Ftry(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #
    ltx    =  ekxx*cotd + ekyx*sotd + cb*lK - cot*Rg
    lty    =  ekxy*cotd + ekyy*sotd         - sot*Rg
    ltz    =  ekxz*cotd             + sb*lK
    #
    return (ltx*ekyx+lty*ekyy+ltz*ekyz)/(ltx*ekzx+lty*ekzy+ltz*ekzz)
    
def Ftrysum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Ftry(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Ftry(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Ftry(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Ftry(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Marx(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #    
    eaxx   =  ekxx*cotd + ekyx*sotd    
    eaxy   =  ekxy*cotd + ekyy*sotd
    eaxz   =  ekxz*cotd
    #
    ltx    =  eaxx + cb*lK - cot*Rg
    lty    =  eaxy         - sot*Rg
    ltz    =  eaxz + sb*lK
    #
    vcx    =  eaxy*ltz - eaxz*lty
    vcy    =  eaxz*ltx - eaxx*ltz
    vcz    =  eaxx*lty - eaxy*ltx
    #
    return (vcx*ekxx+vcy*ekxy+vcz*ekxz)/(ltx*ekzx+lty*ekzy+ltz*ekzz)
    
def Marxsum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Marx(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Marx(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Marx(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Marx(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Mary(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #    
    eaxx   =  ekxx*cotd + ekyx*sotd    
    eaxy   =  ekxy*cotd + ekyy*sotd
    eaxz   =  ekxz*cotd
    #
    ltx    =  eaxx + cb*lK - cot*Rg
    lty    =  eaxy         - sot*Rg
    ltz    =  eaxz + sb*lK
    #
    vcx    =  eaxy*ltz - eaxz*lty
    vcy    =  eaxz*ltx - eaxx*ltz
    vcz    =  eaxx*lty - eaxy*ltx
    #
    return (vcx*ekyx+vcy*ekyy+vcz*ekyz)/(ltx*ekzx+lty*ekzy+ltz*ekzz)
    
def Marysum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Mary(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Mary(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Mary(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Mary(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Marz(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #    
    eaxx   =  ekxx*cotd + ekyx*sotd    
    eaxy   =  ekxy*cotd + ekyy*sotd
    eaxz   =  ekxz*cotd
    #
    ltx    =  eaxx + cb*lK - cot*Rg
    lty    =  eaxy         - sot*Rg
    ltz    =  eaxz + sb*lK
    #
    vcx    =  eaxy*ltz - eaxz*lty
    vcy    =  eaxz*ltx - eaxx*ltz
    vcz    =  eaxx*lty - eaxy*ltx
    #
    return (vcx*ekzx+vcy*ekzy+vcz*ekzz)/(ltx*ekzx+lty*ekzy+ltz*ekzz)
    
def Marzsum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Marz(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Marz(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Marz(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Marz(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Mgrz(lK, Rg, alpha, betas, beta, delta, omegat, ot):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cot    =  np.cos(np.radians(omegat)+ot)
    sot    =  np.sin(np.radians(omegat)+ot)
    cotd   =  np.cos(np.radians(omegat)+ot+delta)
    sotd   =  np.sin(np.radians(omegat)+ot+delta)
    #
    ekxx   =  ca*cbs
    ekxy   =  ca*sbs
    ekxz   = -sa
    #
    ekyx   = -sbs
    ekyy   =  cbs
    ekyz   =  0
    #
    ekzx   =  sa*cbs
    ekzy   =  sa*sbs
    ekzz   =  ca
    #    
    eaxx   =  ekxx*cotd + ekyx*sotd    
    eaxy   =  ekxy*cotd + ekyy*sotd
    eaxz   =  ekxz*cotd
    #
    ltx    =  eaxx + cb*lK - cot*Rg
    lty    =  eaxy         - sot*Rg
    ltz    =  eaxz + sb*lK
    #
    vcz    =  cot*lty - sot*ltx
    #
    return vcz/(ltx*ekzx+lty*ekzy+ltz*ekzz)
    
def Mgrzsum(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    return (Mgrz(lK, Rg, alpha, betas, beta, delta, omegat, pos[0]) +
            Mgrz(lK, Rg, alpha, betas, beta, delta, omegat, pos[1]) +
            Mgrz(lK, Rg, alpha, betas, beta, delta, omegat, pos[2]) +
            Mgrz(lK, Rg, alpha, betas, beta, delta, omegat, pos[3]))/4
    
def Ftrminmax(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    a = np.amin(Ftrxsum(lK, Rg, alpha, betas, beta, delta, omegat, pos))
    b = np.amax(Ftrxsum(lK, Rg, alpha, betas, beta, delta, omegat, pos))
    c = np.amin(Ftrysum(lK, Rg, alpha, betas, beta, delta, omegat, pos))
    d = np.amax(Ftrysum(lK, Rg, alpha, betas, beta, delta, omegat, pos))
    return (np.abs(b-a)+np.abs(d-c)+np.abs(a+b)+np.abs(c+d))
            
def fopt(f, lK, Rg, alpha, betas, beta, delta, omegat, pos):
    def foptfixed(p):
        return f(lK, Rg, p[0], p[1], beta, delta, omegat, pos)
    return foptfixed
    

def orientation(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    # It is very important to start the iterative solution with alpha=0.0 as initial guess 
    # such that the iterative optimizer approaches the solution from the "good" side (horizontal rotor).
    # Especially for configurations at the validity limit with tethers being perpendicular to the
    # rotor axis this is crucial for stable iterative solution.
    guess = [ 0.0, betas ]
    return fmin(fopt(Ftrminmax, lK, Rg, alpha, betas, beta, delta, omegat, pos), guess, xtol=0.00000001, ftol=0.00000001, maxiter=1000)
    
# mean aerodynamic moment over one cycle
def Marzmean(lK, Rg, alpha, betas, beta, delta):
    pos          = np.array([ 0, 0.5*np.pi, np.pi, 1.5*np.pi ])
    omegat       = np.arange(0, 360, 1)
    alpha, betas = orientation(lK, Rg, alpha, betas, beta, delta, omegat, pos)
    return alpha, np.mean(Marzsum(lK, Rg, alpha, betas, beta, delta, omegat, pos)), Ftrminmax(lK, Rg, alpha, betas, beta, delta, omegat, pos)
    
vMarzmean  = np.vectorize(Marzmean)

#f, (ax1, ax2) = plt.subplots(1, 2)

alpha  = np.radians(60)
beta   = np.radians(30)
betas  = np.radians(0)
delta  = np.radians(45)
N      = 4

lK     = np.arange(1.0, 5.2, 0.02) # 0.02 is for print quality
Rg     = np.arange(0.0, 3.2, 0.02) # 0.02 is for print quality
x, y   = np.meshgrid(lK, Rg)

a, z, f   = vMarzmean(x, y, alpha, betas, beta, delta)
h      = np.sin(beta)*lK - np.sin(a)
h1     = np.sin(beta)*lK - np.sin(a)*1.5
h2     = np.sin(beta)*lK - np.sin(a)*2.0
h3     = np.sin(beta)*lK - np.sin(a)*2.5

# apply smoothing - these variables have higher resolution!
zm     = scipy.ndimage.zoom(z, 3)
xx     = scipy.ndimage.zoom(x, 3)
yy     = scipy.ndimage.zoom(y, 3)
ff     = scipy.ndimage.zoom(f, 3)
hh     = scipy.ndimage.zoom(h, 3)
hh1    = scipy.ndimage.zoom(h1, 3)
hh2    = scipy.ndimage.zoom(h2, 3)
hh3    = scipy.ndimage.zoom(h3, 3)

# mask h<0 or too large transverse forces
ffmax  = 10.0
zm     = ma.masked_array(zm, mask=(np.logical_or((hh < 0.0), (ff > ffmax))))
hh     = ma.masked_array(hh, mask=(ff > ffmax))
hh1    = ma.masked_array(hh1, mask=(ff > ffmax))
hh2    = ma.masked_array(hh2, mask=(ff > ffmax))
hh3    = ma.masked_array(hh3, mask=(ff > ffmax))

fig, ax = plt.subplots()

levelsf   = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3]
levels    = [ 0.0,        0.1,         0.2,         0.3,        0.5,         0.7,        1.0]
locations = [(3.0, 0.0), (2.95, 0.5), (2.89, 0.9), (2.8, 1.3), (2.65, 1.8), (2.5, 2.2), (2.4, 2.5)]
loc       = [(1.7, 1.05)]
loc1      = [(2.6, 1.6)]
loc2      = [(3.4, 2.1)]
loc3      = [(4.3, 2.6)]
CS0 = plt.contourf(xx, yy, zm, alpha=0.6, cmap=cm.Blues, levels=levelsf)
CS1 = plt.contour(xx, yy, zm, colors='k', levels=levels)
CS2 = plt.contour(xx, yy, hh, colors='r', linestyles='dotted', levels=[0.0])
CS3 = plt.contour(xx, yy, hh1, colors='r', linestyles='dotted', levels=[0.0])
CS4 = plt.contour(xx, yy, hh2, colors='r', linestyles='dotted', levels=[0.0])
CS5 = plt.contour(xx, yy, hh3, colors='r', linestyles='dotted', levels=[0.0])
CS6 = plt.contour(xx, yy, ff, colors='b', linestyles='dashed', levels=[10.0])
plt.clabel(CS1, inline=1, fontsize=14, manual=locations, fmt='%1.1f')
plt.clabel(CS2, inline=1, fontsize=14, manual=loc, fmt=FormatFaker('1.0'))
plt.clabel(CS3, inline=1, fontsize=14, manual=loc1, fmt=FormatFaker('1.5'))
plt.clabel(CS4, inline=1, fontsize=14, manual=loc2, fmt=FormatFaker('2.0'))
plt.clabel(CS5, inline=1, fontsize=14, manual=loc3, fmt=FormatFaker('2.5'))
#plt.clabel(CS6, inline=1, fontsize=14, manual=loc3, fmt=FormatFaker('limit'))
# reference
plt.plot([1.732,1.732], [0.0,3.0],  color='0.6', linestyle='--')
# high-resolution limit calculation with delta_omega = 0.1
x = [1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60]
y = [1.15, 1.28, 1.38, 1.50, 1.61, 1.73, 1.84, 1.96, 2.07, 2.19, 2.31, 2.54, 2.77, 2.99]
M = [1.52, 1.50, 1.42, 0.97, 1.34, 0.98, 1.27, 1.25, 1.22, 1.21, 1.20, 1.67, 1.14, 1.12]
plt.plot(x, y,  color='k', linestyle='-')
plt.plot([0.0,3.0], [0.0,3.0/np.cos(beta)],  color='r', linestyle='--')
plt.plot([1.732,1.732], [0.0,3.0],  color='0.6', linestyle='--')
plt.scatter(5, 1.25, marker='o', s=80, linewidth=1.2, edgecolor='k', facecolor='w')
plt.xlim(1.0,5.0)
plt.ylim(0.0,3.0)
plt.xticks(np.arange(1.0,5.01,1.0))
plt.yticks(np.arange(0.0,3.01,1.0))
minorLocator = MultipleLocator(0.5)
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_minor_locator(minorLocator)
plt.xlabel(r"Relative distance between rotor centers \$\lK/\Rk\, [-]\$",labelpad=10)
plt.ylabel(r"Rotor size ratio \$\Rg/\Rk\, [-]\$",labelpad=10)

# Test to compare with Fig. 17
a, Ma, f = vMarzmean(5.0, 1.25, alpha, betas, beta, delta)
print "Ma = ", Ma, " alpha = ", np.degrees(a), " betas = ",np.degrees(betas), "Rko/Rk = ",lK*np.sin(beta)/np.sin(a)

# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig18_diagram.svg")
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
