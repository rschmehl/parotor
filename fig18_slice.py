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

dmin = 1.e-5
  
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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, (ltx*ekxx+lty*ekxy+ltz*ekxz)/d, 0.0 )

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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, (ltx*ekyx+lty*ekyy+ltz*ekyz)/d, 0.0 )
    
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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, (vcx*ekxx+vcy*ekxy+vcz*ekxz)/d, 0.0 )
    
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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, (vcx*ekyx+vcy*ekyy+vcz*ekyz)/d, 0.0 )
    
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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, (vcx*ekzx+vcy*ekzy+vcz*ekzz)/d, 0.0 )
    
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
    d      =  ltx*ekzx+lty*ekzy+ltz*ekzz
    #
    return np.where( d > dmin, vcz/d, 0.0 )
    
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
    guess = [ alpha, betas ]
    return fmin(fopt(Ftrminmax, lK, Rg, alpha, betas, beta, delta, omegat, pos), guess, xtol=0.00000001, ftol=0.00000001, maxiter=10000)
    
# mean aerodynamic moment over one cycle
def Marzmean(lK, Rg, alpha, betas, beta, delta):
    pos    = np.array([ 0, 0.5*np.pi, np.pi, 1.5*np.pi ])
    omegat = np.arange(0, 360, 0.1)
    a, b   = orientation(lK, Rg, alpha, betas, beta, delta, omegat, pos)
    alpha  = a
    betas  = b
    Ma     = np.mean(Marzsum(lK, Rg, a, b, beta, delta, omegat, pos))
    Mg     = np.mean(Marzsum(lK, Rg, a, b, beta, delta, omegat, pos))
    Ftau   = Ftrminmax(lK, Rg, alpha, betas, beta, delta, omegat, pos)
    return alpha, betas, Ma, Mg, Ftau
    
vMarzmean  = np.vectorize(Marzmean)

f, (ax1, ax2) = plt.subplots(1, 2)

alpha  = np.radians(60)
betas  = np.radians(0)
beta   = np.radians(30)
delta  = np.radians(45)
#lK     = 1./np.tan(beta)
lK=2.6
Rg     = np.arange(2.6,3.0,0.01)
N      = 4

a, b, Ma, Mg, Ftau = vMarzmean(lK, Rg, alpha, betas, beta, delta)
Mam = np.ma.masked_where(Ftau > 10., Ma)
ax1.plot(Rg, np.degrees(a), 'r')
ax1.plot(Rg, np.degrees(b), 'b')
ax1.plot([0.0,3.0], [60.0,60.0], 'k', label="positive ground clearance")
ax2.plot(Rg, Mam, 'r')
ax2.plot(Rg, Ftau, 'g')
print "lK = ", lK, " Rg = ",Rg[np.argmax(Ftau>10.)-1], " Mam = ",Mam[np.argmax(Ftau>10.)-1]

ax1.set_xlim(0.0,3.0)
ax1.set_ylim(0.0,90.0)
        
ax2.set_xlim(0.0,3.0)
ax2.set_ylim(0.0,2.0)
ax2.yaxis.tick_right()
ax2.yaxis.set_ticks_position('both')

#
# Set axis labels
plt.title(r"\$\lK/\Rk\=$"+str(lK))
ax1.set_xlabel(r"Rotor size ratio \$\Rg/\Rk, [-]\$",labelpad=10)
ax1.set_ylabel(r"Flow angles \$[-]\$",labelpad=10)
ax1.legend(loc='upper right', frameon=False)
ax2.set_xlabel(r"Rotor size ratio \$\Rg/\Rk, [-]\$",labelpad=10)
ax2.set_ylabel(r"Transferred moment \$[-]\$",labelpad=10)
ax2.yaxis.set_label_position("right")
ax2.legend(loc='bottom center', frameon=False)
#
# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig18_slice.svg")
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    