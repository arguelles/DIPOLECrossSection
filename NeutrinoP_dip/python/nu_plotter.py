from numpy import sqrt, genfromtxt
from pylab import subplots, cla, clf, close, savefig
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, LogLocator
from scipy.interpolate import interp1d
import re

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'], 'size' : 18})
rc('legend', fontsize = 14)

#       blue     | orange  | green   | yellow  | magenta | grey
cols = ['#29A2C6','#FF6D31','#73B66B','#FFCB18','#EF597B', '#333333']
#       blue     | orange  | yellow  | green   | magenta | grey
#cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']
#       blue     | green   | magenta | grey
#cols = ['#29A2C6','#73B66B','#EF597B', '#333333']

#BBOX-Properties
bbox_props = dict(boxstyle="round,pad=0.3", fc='#FFF6EA', ec="black", lw=0.3)
bbox_props1 = dict(boxstyle="round,pad=0.3", fc='#FFFFFF', ec="white", lw=0.5)

root="/home/mkroll/work/GammaP_dip/"
root_table= root + "data/"

dipole_display = {
		'gbw'       : 'GBW (1999)',
		'machado'   : 'Machado (2008)',
		'maddip'    : r'\textbf{MadDip}',
		'munier'    : 'Munier',
		'soyez'     : 'Soyez (2007)',
}

dipole_colors = {
		'gbw'       : cols[2],
		'machado'   : cols[5],
		'maddip'    : cols[1],
		'munier'    : cols[0],
		'soyez'     : cols[4],
}


dipole_style= {
		'gbw'       : '-.',
		'machado'   : 'solid',
		'maddip'    : 'solid',
		'munier'    : '--',
		'soyez'     : ':',
}

ext_display = {
        'cc'        : ' CC',
        'nc'        : ' NC',
}

#def Plot(reference, fits):
def Plot(fits, extension):
    print "Plot parameters "
    fig, ax = subplots(nrows=1, ncols=1, squeeze=True, 
        sharex=False, sharey=False, figsize=(5*(1.+sqrt(5.))/2., 5))
    ax.set_xlabel(r'$E_{\gamma}$ [GeV]')
    ax.set_ylabel(r'$\sigma_{\gamma P}$ [cm$^2$]')
    ax.grid(True)
    i=0
   # data_ref = genfromtxt(reference)
   # x_ref = data_ref[:,0]
   # y_ref = data_ref[:,1]

   # ax.loglog(x_ref, y_ref, color = 'k', lw = 2, 
   #     label = "Machado (2008)", alpha = 0.6)

    for fit in fits:
        for ext in extension: 
            file = root_table + "neutrino/" + fit + "_" + ext + ".dat" 
            data = genfromtxt(file, skip_header=3)
            x_val = data[:,0]
            y_val = data[:,1]
            if (ext == 'nc'):
                ax.loglog(x_val, y_val, color = dipole_colors[fit], ls = ':', lw = 3, 
                    label = dipole_display[fit] + ext_display[ext], alpha = 1.0)
            else:
                ax.loglog(x_val, y_val, color = dipole_colors[fit], ls = dipole_style[fit], lw = 3, 
                    label = dipole_display[fit] + ext_display[ext], alpha = 1.0)
                #label = pdf_display[pdfs] + ' (CC)', alpha = 1.0)
            i += 1
            x_val = []
            y_val = []

    i=0

    ax.legend(loc = 4, fancybox = True, shadow = True, numpoints = 1)
    ax.set_ylim(1.e-37, 1.e-31)
    ax.set_xlim(1.e3, 2.e12)
    savefig("/home/mkroll/work/GammaP_dip/plots/neutrino.pdf", bbox_inches='tight')
    cla()
    clf()
    close()


if __name__ == "__main__":
    data_list = []
    dipoles = ["maddip"]
    ext = ["cc", "nc"]
    #pdfsets = ["CT10nlo", "HERAPDF15NLO_EIG"]
    reference = '/home/mkroll/work/GammaP_dip/data/machado.dat'

    #data_list.append(genfromtxt('/home/mkroll/nu_XS/totalXS_Sakar.dat'))

    #for pdf in pdfsets:
    #    data_list.append(genfromtxt(root_table + pdf + '_nu_CCxs.dat', skip_header=2))

    #Plot(reference, dipoles)
    Plot(dipoles, ext)
