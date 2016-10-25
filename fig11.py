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

# tether length without suspension lines (alpha and beta independent)
def ltfree(lK, Rg, alpha, beta, delta, omegat):
    ltx    =  np.cos(alpha)*np.cos(np.radians(omegat)+delta) + np.cos(beta)*lK - np.cos(np.radians(omegat))*Rg
    lty    =                np.sin(np.radians(omegat)+delta)                   - np.sin(np.radians(omegat))*Rg
    ltz    = -np.sin(alpha)*np.cos(np.radians(omegat)+delta) + np.sin(beta)*lK
    return np.sqrt(ltx*ltx+lty*lty+ltz*ltz)

# min tether length over cycle
def ltfreemin(lK, Rg, alpha, beta, delta): 
    omegat = np.arange(0, 360, 1)
    return np.amin(ltfree(lK, Rg, alpha, beta, delta, omegat))

# max tether length over cycle
def ltfreemax(lK, Rg, alpha, beta, delta): 
    omegat = np.arange(0, 360, 1)
    return np.amax(ltfree(lK, Rg, alpha, beta, delta, omegat))
    
# tether length
def lt(lK, Rg, beta, delta, omegat):
    alpha  = 0.5*np.pi-beta
    ltx    =  np.cos(alpha)*np.cos(np.radians(omegat)+delta) + np.cos(beta)*lK - np.cos(np.radians(omegat))*Rg
    lty    =                np.sin(np.radians(omegat)+delta)                   - np.sin(np.radians(omegat))*Rg
    ltz    = -np.sin(alpha)*np.cos(np.radians(omegat)+delta) + np.sin(beta)*lK
    return np.sqrt(ltx*ltx+lty*lty+ltz*ltz)
    
# min tether length over cycle
def ltmin(lK, Rg, beta, delta): 
    omegat = np.arange(0, 360, 1)
    return np.amin(lt(lK, Rg, beta, delta, omegat))

# max tether length over cycle
def ltmax(lK, Rg, beta, delta): 
    omegat = np.arange(0, 360, 1)
    return np.amax(lt(lK, Rg, beta, delta, omegat))

vltfreemin = np.vectorize(ltfreemin)
vltfreemax = np.vectorize(ltfreemax)
vltmin = np.vectorize(ltmin)
vltmax = np.vectorize(ltmax)

f, (ax1, ax2) = plt.subplots(1, 2)

lK    = 2
Rg    = 1
beta  = np.radians(np.arange(0, 91, 1) )
delta = np.radians(0)
ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
#print 'ltmin & ltmax ', ltmin(lK, Rg, np.radians(30), delta), ltmax(lK, Rg, np.radians(30), delta)
delta = np.radians(30)
ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
delta = np.radians(60)
ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
delta = np.radians(90)
ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
delta = np.radians(120)
#ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
#ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
delta = np.radians(150)
#ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
#ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
delta = np.radians(180)
#ax1.plot(np.degrees(beta), vltmax(lK, Rg, beta, delta), 'r')
#ax1.plot(np.degrees(beta), vltmin(lK, Rg, beta, delta), 'b')
plt.xlim(0,90)
ax1.set_xticks(np.arange(0,91,15))
ax1.set_ylim(0.5,3.5)

lK    = 2
Rg    = 1
beta  = np.radians(30)
alpha = np.radians(np.arange(0, 91, 1) )
delta = np.radians(0)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
#print 'ltmin & ltmax ', ltfreemin(lK, Rg, np.radians(60), beta, delta), ltfreemax(lK, Rg, np.radians(60), beta, delta)
delta = np.radians(30)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
delta = np.radians(60)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
delta = np.radians(90)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
delta = np.radians(120)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
delta = np.radians(150)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
delta = np.radians(180)
ax2.plot(np.degrees(alpha), vltfreemax(lK, Rg, alpha, beta, delta), 'r')
ax2.plot(np.degrees(alpha), vltfreemin(lK, Rg, alpha, beta, delta), 'b')
ax2.set_xticks(np.arange(0,91,15))
ax2.set_ylim(0.5,3.5)

#
# Set axis labels
ax1.set_xlabel(r"Elevation angle \$\beta\, [^{\circ}]\$",labelpad=10)
ax2.set_xlabel(r"Inclination angle \$\alpha\, [^{\circ}]\$",labelpad=10)
ax1.set_ylabel(r"Tether length \$\lt/\Rk\, [-]\$",labelpad=10)

#
# Save as SVG file
plt.subplots_adjust(left=0.07, bottom=None, right=0.93, top=None, wspace=None, hspace=None)
plt.savefig("fig11_diagram.svg")
plt.show()
