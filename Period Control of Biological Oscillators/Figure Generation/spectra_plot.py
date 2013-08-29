# Produces Figure 5

import os
from analysis import *

Plotting.figure(figsize=(6,6))

# Load the eigvals
c = 10
e=list()
e.append(control_direction('SubstrateDepletion_inc_ens_final')[1][0])
e.append(control_direction('SubstrateDepletion_dec_ens_final')[1][0])
e.append(control_direction('repressilator_inc_ens_final')[1][0])
e.append(control_direction('repressilator_dec_ens_final')[1][0])
e.append(control_direction('repressilator_positive_inc_ens_final')[1][0])
e.append(control_direction('repressilator_positive_dec_ens_final')[1][0])
e.append(control_direction('goodwin_inc_ens_final')[1][0])
e.append(control_direction('goodwin_dec_ens_final')[1][0])
e.append(control_direction('calcium_oscillations_inc_ens_final')[1][0])
e.append(control_direction('calcium_oscillations_dec_ens_final')[1][0]) #9

for model_idx, eigen in enumerate(e):
    width = (234/len(eigen))**0.25 * 0.25
    l = Plotting.plot_eigval_spectrum(eigen/max(eigen), offset = 0.15+model_idx, 
                                      widths=0.7, lw=width)

# Now a lot of fiddling to make the plot prettier
ax = Plotting.gca()
for ii in range(1, c):
    ax.plot([ii, ii], [0.5e-6, 2], '-', lw=0.5, color=(0.75,0.75,0.75))

# Add labels
ax.set_xticks(0.5 + scipy.arange(c))
import string
xlabels = [r'$SD_{inc}$', r'$SD_{dec}$', r'$R_{inc}$', r'$R_{dec}$',
           r'$RP_{inc}$', r'$RP_{dec}$', r'$G_{inc}$', r'$G_{dec}$', 
           r'$C_{inc}$',r'$C_{dec}$']

#['%s' % tup for tup in string.ascii_lowercase[:c]]                                          
ax.set_xticklabels(xlabels, fontsize=12, verticalalignment='bottom',
                   rotation=0, horizontalalignment='center')
for l in ax.get_xticklabels():
    l.set_y(l.get_position()[1] - 0.06)

for l in ax.get_xticklines():
    l.set_visible(False)

ticks = range(-6, 1)
ax.set_yticks([10**ii for ii in ticks])
import matplotlib.lines

for l in ax.get_yticklines():
    if l.get_xdata() == (1,):
        l.set_visible(False)
    l.set_marker(matplotlib.lines.TICKLEFT)

ax.set_yticklabels([r'$10^{%s}$' % ii for ii in ticks], fontsize=12)
ax.set_ylabel(r'$\lambda / \lambda_{max}$', fontsize=12)
ax.set_title('PCA Eigenvalue Hierarchy')


#ax.set_xlim(0, c)

#ax.set_yscale('log',subsy=[])
#for l in ax.get_yticklabels():
#   l.set_x(l.get_position()[0] - 0.04)

ax.set_ylim(0.5e-6, 2)
#Plotting.subplots_adjust(bottom=0.08, right=0.97, top=0.97)

Plotting.savefig('eig_spectra_plot.eps')
