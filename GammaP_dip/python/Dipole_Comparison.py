#!/usr/bin/python

from pylab import *
from matplotlib import rc, rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy.interpolate import interp1d

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'], "size": 14})
#rc('legend', fontsize = 16)

# Read in data from an ASCII data table
#data103 = genfromtxt('../data/dipole/dipole_comparison.dat')
data103 = genfromtxt('/home/mkroll/work/GammaP_dip/data/dipole/dipole_comparison.dat')

#'data' is a matrix containing the columns and rows from the file
tau     = data103[:,0]  # Python indices are (row,col) as in linalg
sig_ren = data103[:,1]  # Creates arrays for first two columns
sig_cgc = data103[:,2]
sig_gbw = data103[:,3]
sig_soy = data103[:,4]


cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B']

fig, ax = subplots(nrows=1, ncols=1, squeeze=True, sharex=False, sharey=False, figsize=((1.+sqrt(5.))/2.*5, 5))

# Create a plot of data
ax.plot(tau, sig_ren, color = cols[1], ls = 'solid', lw = 3, label = r'\textbf{This work}')
ax.plot(tau, sig_cgc, color = cols[0], ls = '--', lw = 2, label = r'Munier (2004)')
ax.plot(tau, sig_gbw, color = cols[3], ls = '-.', lw = 2, label = r'GBW (1999)')
ax.plot(tau, sig_soy, color = cols[4], ls = ':', lw = 2, label = r'Soyez (2007)')

# Set label
ax.set_xlabel(r'$\tau = r Q_s = r Q_0 \left( \frac{x_0}{x} \right)^{\lambda/2}$')
ax.set_ylabel(r'$\sigma (\tau)$ [mb]')

# set legend
ax.legend(loc = 2, fancybox = True, ncol = 2, shadow = True)

ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
ax.yaxis.set_minor_locator(AutoMinorLocator())

bbox_props = dict(boxstyle="round,pad=0.3", fc='#FFF6EA', ec="black", lw=0.3)
ax.text(8.4, 2.5, r'$x = 10^{-5}$', bbox=bbox_props, fontsize = 16)

# Turn on a grid
#grid(True)

#fig.tight_layout()

# Set axes range
ax.set_ylim([0,50])

# Set plot title
#title(r'Comparison of Dipole Models at $x = 10^{-5}$')

# Save the figure in a separate file
fig.savefig('/home/carguelles/work/GammaP_dip/plots/Dipole_Comparison.pdf', bbox_inches='tight')

# Draw the plot to the screen
#fig.show()
#raw_input()
