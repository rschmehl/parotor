from mpl_toolkits.mplot3d import axes3d, proj3d
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import np
###patch start### http://stackoverflow.com/questions/16488182/removing-axes-margins-in-3d-plot
from mpl_toolkits.mplot3d.axis3d import Axis
if not hasattr(Axis, "_get_coord_info_old"):
    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs
    Axis._get_coord_info_old = Axis._get_coord_info  
    Axis._get_coord_info = _get_coord_info_new
###patch end###

# Figure size
mpl.rcParams['figure.figsize'] = 10, 10.5

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

def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return np.array([[1,0,0,0],
                     [0,1,0,0],
                     [0,0,a,b],
                     [0,0,-0.0001,zback]])

# Propblem parameters
beta   = np.radians(30)
alpha  = np.radians(60) # 90 - beta
betas  = np.radians(10)
delta  = np.radians(30)
lK     = 100
Rg     = 25
Rk     = 20
Rko    = 35
le     = 40
N      = 4
ot0    = 40 # offset
ot     = np.array([ 0, 90, 180, 270 ])

# Derived trigonometric data
omegat = np.arange(0, 361, 1)

def rB(R, omegat):
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    rBx    = cot*R
    rBy    = sot*R
    rBz    = 0
    return np.array([ rBx, rBy, rBz ])
    
def rA(lK, R, alpha, betas, beta, delta, omegat):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    rAx    =  (ca*cbs*cotd - sbs*sotd)*R + cb*lK
    rAy    =  (ca*sbs*cotd + cbs*sotd)*R
    rAz    =  (   -sa*cotd           )*R + sb*lK
    return np.array([ rAx, rAy, rAz ])
    
def rK(lK, beta):
    cb     =  np.cos(beta)
    sb     =  np.sin(beta)
    rKx    =  cb*lK
    rKy    =  0
    rKz    =  sb*lK
    return np.array([ rKx, rKy, rKz ])

def ew():
    return np.array([[ 1, 0, 0], 
                     [ 0, 1, 0], 
                     [ 0, 0, 1]])
                     
def ek(alpha, betas):
    ca     =  np.cos(alpha)
    sa     =  np.sin(alpha)
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)
    return np.array([[ ca*cbs, -sbs,  sa*cbs], 
                     [ ca*sbs,  cbs,  sa*sbs], 
                     [-sa,      0,    ca]])
                     
def es(betas):
    cbs    =  np.cos(betas)
    sbs    =  np.sin(betas)   
    return np.array([[ cbs, -sbs,  0], 
                     [ sbs,  cbs,  0], 
                     [   0,    0,  1]])
                     
def ea(alpha, betas, delta, omegat):
    cotd   =  np.cos(np.radians(omegat)+delta)
    sotd   =  np.sin(np.radians(omegat)+delta)
    T      =  ek(alpha, betas)
    e      =  np.array([[ cotd,-sotd, 0],
                        [ sotd, cotd, 0],
                        [    0,    0, 1]])
    return T.dot(e)

def eb(omegat):
    cot    =  np.cos(np.radians(omegat))
    sot    =  np.sin(np.radians(omegat))
    return np.array([[ cot,-sot, 0],
                     [ sot, cot, 0],
                     [   0,   0, 1]])
                     
def arc(r0, R, e, n, phi0, phi):
    e1 = e/np.sqrt(np.sum(e*e)) # normalize
    en = n/np.sqrt(np.sum(n*n)) # normalize
    ip = np.argmax(phi>phi0)    # find end index
    e2 = np.cross(en, e1)
    cp = np.cos(np.radians(phi[:ip]))
    sp = np.sin(np.radians(phi[:ip]))
#    r  = cp*e1+sp*e2
    r = np.zeros((3, ip))
    r[0,:] = r0[0]+R*(cp*e1[0]+sp*e2[0])
    r[1,:] = r0[1]+R*(cp*e1[1]+sp*e2[1])
    r[2,:] = r0[2]+R*(cp*e1[2]+sp*e2[2])
    return r
    
fig = plt.figure()
ax  = fig.gca(projection='3d')
#ax.set_axis_off() 
#ax.w_zaxis.line.set_lw(0.)
#ax.set_zticks([])
for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)
ax.w_zaxis.gridlines.set_visible(False)
ax.w_zaxis.pane.set_visible(False)
ax.w_zaxis.line.set_visible(False)
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

# Ground rotor
r = rB(Rg, omegat)
ax.plot(r[0], r[1], r[2], linewidth=1.2, color='b')
ax.scatter(0, 0, 0, zorder=20, s=100, linewidth=1.2, edgecolor='b', facecolor='b')
r = rB(Rg, ot0+ot[0])
ax.scatter(r[0], r[1], r[2], zorder=20, s=100, linewidth=1.2, edgecolor='b', facecolor='w')
r = rB(Rg, ot0+ot[1])
ax.scatter(r[0], r[1], r[2], zorder=20, s=100, linewidth=1.2, edgecolor='b', facecolor='w')
r = rB(Rg, ot0+ot[2])
ax.scatter(r[0], r[1], r[2], zorder=20, s=100, linewidth=1.2, edgecolor='b', facecolor='w')
r = rB(Rg, ot0+ot[3])
ax.scatter(r[0], r[1], r[2], zorder=20, s=100, linewidth=1.2, edgecolor='b', facecolor='w')

# vector base wind reference frame
e = le*ew()
ewx = e[:,0]
ax.plot([0, ewx[0]], [0, ewx[1]], [0, ewx[2]], linewidth=3, color='k')
ewy = e[:,1]
ax.plot([0, ewy[0]], [0, ewy[1]], [0, ewy[2]], linewidth=3, color='k')
ewz = e[:,2]
ax.plot([0, ewz[0]], [0, ewz[1]], [0, ewz[2]], linewidth=3, color='k')

# vector base rotating reference frame
e = le*eb(ot0)
ebx = e[:,0]
ax.plot([0, ebx[0]], [0, ebx[1]], [0, ebx[2]], linewidth=3, color='b')
eby = e[:,1]
ax.plot([0, eby[0]], [0, eby[1]], [0, eby[2]], linewidth=3, color='b')
ebz = e[:,2]
#ax.plot([0, ebz[0]], [0, ebz[1]], [0, ebz[2]], linewidth=2, color='b')

# elevation angle
r0   = np.array([0, 0, 0])
R    = 32
e    = ewx
n    =-ewy
r    = arc(r0, R, e, n, np.degrees(beta), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], zorder=20, linewidth=2, color='r')

# rotor spinning
R    = 20
e    = ewx
n    = ewz
r    = arc(r0, R, e, n, ot0, omegat)
ax.plot(r[0,:], r[1,:], r[2,:], zorder=20, linewidth=2, color='r')

# Flying rotor
r = rA(lK, Rk, alpha, betas, beta, delta, omegat)
ax.plot(r[0], r[1], r[2], zorder=-10, linewidth=1.2, color='#00ff00')
r = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[0])
ax.scatter(r[0], r[1], r[2], zorder=10, s=100, linewidth=1.2, edgecolor='#00ff00', facecolor='w')
r = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[1])
ax.scatter(r[0], r[1], r[2], zorder=10, s=100, linewidth=1.2, edgecolor='#00ff00', facecolor='w')
r = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[2])
ax.scatter(r[0], r[1], r[2], zorder=10, s=100, linewidth=1.2, edgecolor='#00ff00', facecolor='w')
r = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[3])
ax.scatter(r[0], r[1], r[2], zorder=10, s=100, linewidth=1.2, edgecolor='#00ff00', facecolor='w')
r = rA(lK, Rko, alpha, betas, beta, delta, omegat)
ax.plot(r[0], r[1], r[2], linewidth=1.2, color='#00ff00')
r = rK(lK, beta)
ax.scatter(r[0], r[1], r[2], zorder=5, s=100, linewidth=1.2, edgecolor='#00ff00', facecolor='#00ff00')

# vector base wind reference frame at kite location
ax.plot([r[0], r[0]+ewx[0]], [r[1], r[1]+ewx[1]], [r[2], r[2]+ewx[2]], zorder=-20, linewidth=3, color='0.7')
ax.plot([r[0], r[0]+ewy[0]], [r[1], r[1]+ewy[1]], [r[2], r[2]+ewy[2]], zorder=-20, linewidth=3, color='0.7')
#ax.plot([r[0], r[0]+ewz[0]], [r[1], r[1]+ewz[1]], [r[2], r[2]+ewz[2]], zorder=-20, linewidth=2, color='0.7')

# vector base intermediate reference frame
e = le*es(betas)
esx = e[:,0]
ax.plot([r[0], r[0]+esx[0]], [r[1], r[1]+esx[1]], [r[2], r[2]+esx[2]], linewidth=3, color='0.7')
esy = e[:,1]
ax.plot([r[0], r[0]+esy[0]], [r[1], r[1]+esy[1]], [r[2], r[2]+esy[2]], linewidth=3, color='0.7')
esz = e[:,2]
ax.plot([r[0], r[0]+esz[0]], [r[1], r[1]+esz[1]], [r[2], r[2]+esz[2]], linewidth=3, color='0.7')

# vector base kite reference frame
e = le*ek(alpha, betas)
ekx = e[:,0]
ax.plot([r[0], r[0]+ekx[0]], [r[1], r[1]+ekx[1]], [r[2], r[2]+ekx[2]], linewidth=3, color='k')
eky = e[:,1]
ax.plot([r[0], r[0]+eky[0]], [r[1], r[1]+eky[1]], [r[2], r[2]+eky[2]], linewidth=3, color='k')
ekz = e[:,2]
ax.plot([r[0], r[0]+ekz[0]], [r[1], r[1]+ekz[1]], [r[2], r[2]+ekz[2]], linewidth=3, color='k')

# vector base kite rotating reference frame
e = le*ea(alpha, betas, delta, ot0)
eax = e[:,0]
ax.plot([r[0], r[0]+eax[0]], [r[1], r[1]+eax[1]], [r[2], r[2]+eax[2]], zorder=-20, linewidth=3, color='#00ff00')
eay = e[:,1]
ax.plot([r[0], r[0]+eay[0]], [r[1], r[1]+eay[1]], [r[2], r[2]+eay[2]], zorder=-20, linewidth=3, color='#00ff00')
eaz = e[:,2]
#ax.plot([r[0], r[0]+eaz[0]], [r[1], r[1]+eaz[1]], [r[2], r[2]+eaz[2]], zorder=-20, linewidth=2, color='#00ff00')

# additional base vectors
#ax.plot([r[0], r[0]+ewz[0]], [r[1], r[1]+ewz[1]], [r[2], r[2]+ewz[2]], linewidth=2, color='0.7')

# lK
ax.plot([0, r[0]], [0, r[1]], [0, r[2]], linewidth=1.22, color='k')

# kite reference frame rotation
r0   = r
R    = 35
e    = ewx
n    = ewz
r    = arc(r0, R, e, n, np.degrees(betas), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
e    = ewy
r    = arc(r0, R, e, n, np.degrees(betas), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')

R    = 35
e    = esx
n    = esy
r    = arc(r0, R, e, n, np.degrees(alpha), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
e    = esz
r    = arc(r0, R, e, n, np.degrees(alpha), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')

# rotor spinning
R    = 15
e    = ekx
n    = ekz
r    = arc(r0, R, e, n, ot0, omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
r    = arc(r0, R, e, n, ot0+np.degrees(delta), omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')

# tethers
ra = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[0])
rb = rB(Rg, ot0+ot[0])
ax.plot([rb[0], ra[0]], [rb[1], ra[1]], [rb[2], ra[2]], zorder=-20, linewidth=1.2, color='k')
#
ra = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[1])
rb = rB(Rg, ot0+ot[1])
ax.plot([rb[0], ra[0]], [rb[1], ra[1]], [rb[2], ra[2]], zorder=-20, linewidth=1.2, color='0.8')
#
ra = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[2])
rb = rB(Rg, ot0+ot[2])
ax.plot([rb[0], ra[0]], [rb[1], ra[1]], [rb[2], ra[2]], zorder=-20, linewidth=1.2, color='0.8')
#
ra = rA(lK, Rk, alpha, betas, beta, delta, ot0+ot[3])
rb = rB(Rg, ot0+ot[3])
ax.plot([rb[0], ra[0]], [rb[1], ra[1]], [rb[2], ra[2]], zorder=-20, linewidth=1.2, color='0.8')

# kite angular velocities
rk0  = rK(lK, beta)
r0   = rk0+0.8*ekz
R    = 5
e    = 0.5*ekx-eky
n    = ekz
r    = arc(r0, R, e, n, 320, omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
ax.scatter(r0[0], r0[1], r0[2], s=10, linewidth=1.2, edgecolor='r', facecolor='w')
r0   = rk0+1.2*ekx
R    = 5
e    =-eky-ekz
n    = ekx
r    = arc(r0, R, e, n, 320, omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
ax.scatter(r0[0], r0[1], r0[2], s=10, linewidth=1.2, edgecolor='r', facecolor='w')
r0   = rk0+1.2*eky
R    = 5
e    =-0.5*ekx+ekz
n    = eky
r    = arc(r0, R, e, n, 320, omegat)
ax.plot(r[0,:], r[1,:], r[2,:], linewidth=2, color='r')
ax.scatter(r0[0], r0[1], r0[2], s=10, linewidth=1.2, edgecolor='r', facecolor='w')



ax.set_xlim3d(-25, 75)
ax.set_ylim3d(-50, 50)
ax.set_zlim3d(  0, 100)
#ax.axis('equal')

#ax.set_aspect('equal')

ax.view_init(20, 225)
ax.set_axis_off()

# orthogonal projection
proj3d.persp_transformation = orthogonal_proj

plt.tight_layout()
plt.subplots_adjust(left=-0.31, right=1.12, top=1.4, bottom=-0.3)
#plt.subplots_adjust(left=-0.31, right=1.12, top=1.2, bottom=-0.5) # to get the top arc
plt.savefig("fig09_diagram.svg")
plt.show()