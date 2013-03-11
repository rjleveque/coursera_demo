
from matplotlib import pyplot as plt
from matplotlib import cm               # colormaps
import numpy as np

try:
    fname = 'solution.txt'
    x, y, u = np.loadtxt(fname, unpack=True)
except:
    err_msg = "Could not load data from file %s." % fname \
              + " Did you forget to run the program?"
    raise Exception(err_msg)


# Solution is plotted on n by n grid so length of each vector should be n**2
# Determine n:

n = int(np.sqrt(len(x)))
assert n*n == len(x), "Expected len(x) to be a perfect square, len(x) = %s" % len(x)


X = x.reshape(n,n)
Y = y.reshape(n,n)
U = u.reshape(n,n)

# Pseudocolor plot

plt.figure(1)                  
plt.clf()                      # clear figure
plt.axis('scaled')             # so x- and y-axis scaled the same (square)
plt.pcolor(X,Y,U,cmap=cm.jet)  # pseudo-color plot using colormap "jet"
plt.clim(0., 1.)               # colors range from u=0 to u=1
plt.colorbar()                 # add a color bar to show temperature scale
plt.title('Temperature')

plt.savefig('pcolor.png')
print 'Saved pseudocolor plot as pcolor.png'

# Contour plot

plt.figure(2)                  
plt.clf()                     
plt.axis('scaled')             

# contour line levels:
clines = np.linspace(0., 1., 26)

# do contour plot:
C = plt.contour(X,Y,U,clines,colors='k') 

# add labels on every other line:
plt.clabel(C, clines[1::2], inline=1, fontsize=10)

plt.title('Contours of temperature')

plt.savefig('contour.png')
print 'Saved contour plot as contour.png'
