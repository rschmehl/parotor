from pylab import np

Ftmin = 0.88
Ftmax = 1.76
#ltmin = 0.6
#ltmax = 1.2
ltmin = 0.8
ltmax = 1.4
omega = 3.0

print "Preel = ", (Ftmax+Ftmin)*(ltmax-ltmin)*omega/np.pi

