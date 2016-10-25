import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import np
from scipy.optimize import fmin
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
#    print np.abs(a+b)+np.abs(c+d)
    return (np.abs(a+b)+np.abs(c+d))
            
def fopt(f, lK, Rg, alpha, betas, beta, delta, omegat, pos):
    def foptfixed(p):
        return f(lK, Rg, p[0], p[1], beta, delta, omegat, pos)
    return foptfixed
    

def orientation(lK, Rg, alpha, betas, beta, delta, omegat, pos):
    guess = [ alpha, betas ] # can be improved by passing proving
    return fmin(fopt(Ftrminmax, lK, Rg, alpha, betas, beta, delta, omegat, pos), guess, xtol=0.00000001, ftol=0.00000001, maxiter=1000)
    
# mean aerodynamic moment over one cycle
def Marzmean(lK, Rg, alpha, betas, beta, delta):
    pos    = np.array([ 0, 0.5*np.pi, np.pi, 1.5*np.pi ])
    omegat = np.arange(0, 360, 1)
    a, b   = orientation(lK, Rg, alpha, betas, beta, delta, omegat, pos)
    alpha  = a
    betas  = b
    Ma     = np.mean(Marzsum(lK, Rg, a, b, beta, delta, omegat, pos))
    Mg     = np.mean(Marzsum(lK, Rg, a, b, beta, delta, omegat, pos))
    return alpha, betas, Ma, Mg
    
def Marzmean0(lK, Rg, alpha, betas, beta, delta):
    pos    = np.array([ 0, 0.5*np.pi, np.pi, 1.5*np.pi ])
    omegat = np.arange(0, 360, 1)
    return np.mean(Marzsum(lK, Rg, alpha, betas, beta, delta, omegat, pos))
    
vMarzmean  = np.vectorize(Marzmean)
vMarzmean0 = np.vectorize(Marzmean0)

f, (ax1, ax2) = plt.subplots(1, 2)

alpha  = np.radians(60) # 90 - beta | 62.3608549873
betas  = np.radians(0)  # 4.38381787735
lK     = 2
Rg     = 1
N      = 4

delta  = np.arange(0, 181, 1)

beta   = np.radians(30)
a, b, Ma, Mg = vMarzmean(lK, Rg, alpha, betas, beta, np.radians(delta))
ax1.plot(delta, np.degrees(a), 'r', label=r"\$\alpha\$")
ax1.plot(delta, np.degrees(b), 'b', label=r"\$\betas\$")
ax2.plot(delta, Ma, 'r', label=r"\$\Mamean/(\Rk\Fakz)\$")
#ax2.plot(delta, vMarzmean0(lK, Rg, alpha, betas, beta, np.radians(delta)), 'b--')  

beta   = np.radians(45)
a, b, Ma, Mg = vMarzmean(lK, Rg, alpha, betas, beta, np.radians(delta))
ax1.plot(delta, np.degrees(a), 'r')
ax1.plot(delta, np.degrees(b), 'b--')
ax2.plot(delta, Ma, 'r')
#ax2.plot(delta, vMarzmean0(lK, Rg, alpha, betas, beta, np.radians(delta)), 'b--')

beta   = np.radians(60)
a, b, Ma, Mg = vMarzmean(lK, Rg, alpha, betas, beta, np.radians(delta))
ax1.plot(delta, np.degrees(a), 'r')
ax1.plot(delta, np.degrees(b), 'b--')
ax2.plot(delta, Ma, 'r')
#ax2.plot(delta, vMarzmean0(lK, Rg, alpha, betas, beta, np.radians(delta)), 'b--')

ax1.set_xlim(0,180)
ax1.set_ylim(0,70)
ax1.set_xticks(np.arange(0,181,30))


#alpha  = np.radians(60)
#betas  = np.radians(0)
#ax2.plot(delta, Ma, 'r')
#ax2.plot(delta, vMarzmean0(lK, Rg, alpha, betas, beta, np.radians(delta)), 'b--')          
ax2.set_xlim(0,180)
ax2.set_ylim(0,0.5)
ax2.set_xticks(np.arange(0,181,30))
ax2.set_yticks(np.arange(0,0.51,0.1))
ax2.yaxis.tick_right()
ax2.yaxis.set_ticks_position('both')

#
# Set axis labels
ax1.set_xlabel(r"Phase lag angle \$\delta\, [^{\circ}]\$",labelpad=10)
ax1.set_ylabel(r"Flow angles \$[-]\$",labelpad=10)
ax1.legend(loc='upper right', frameon=False)
ax2.set_xlabel(r"Phase lag angle \$\delta\, [^{\circ}]\$",labelpad=10)
ax2.set_ylabel(r"Transferred moment \$[-]\$",labelpad=10)
ax2.yaxis.set_label_position("right")
ax2.legend(loc='bottom center', frameon=False)
#
# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig17_diagram.svg")
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    