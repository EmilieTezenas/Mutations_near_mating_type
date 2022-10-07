import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import math
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


#===================
#Computation of the eigenvalues
#===================

def alpha_overdom(r, pin):
    return pin*(1-r)+r-2

def beta_overdom(r, pin, s3, s4, f):

    alpha = alpha_overdom(r, pin)
    return f*(r*alpha - 2*s3 +2) - 2 + 4*s4

def Delta_overdom(r, pin, s3, s4, f):

    alpha = alpha_overdom(r, pin)

    return f**2*(alpha*r - 2*s3 +2)**2 + 4*f*(2*s3 - 2 + r*(2*s3-1)*alpha) + 4

def valeurs_propres_overdom(r, pin, s3, s4, f):

    beta = beta_overdom(r, pin, s3, s4, f)
    Delta = Delta_overdom(r, pin, s3, s4, f)
    racine_delta = np.sqrt(Delta)

    return s4-r, 0.25*(beta + racine_delta), 0.25*(beta - racine_delta)

#===================
#To compute the matrix of size len(Echelle_R)xlen(Echelle_s4) with the value of the dominant eigenvalue for corresponding values of r and s4, for fixed values of f, s3 and pin.
#===================

def evol_vap_overdom_s4_r_nograph(Echelle_R, pin, Echelle_s4, s3, f):

    n = len(Echelle_R)
    m = len(Echelle_s4)
    les_vap = np.zeros((n,m))

    for k in range(n):
        r = Echelle_R[n-1-k]#Pour que l'axe des ordonnées soit dans le sens croissant en "montant"
        for l in range(m):
            s4 = Echelle_s4[l]
            lambda_plus = valeurs_propres_overdom(r, pin, s3, s4, f)[1]
            les_vap[k][l] = lambda_plus

    min = np.min(les_vap)
    max = np.max(les_vap)

    #On met en forme les données pour que seaborn prenne les bons arguments en compte et que ce soit joli
    data0 = pd.DataFrame(les_vap[:][:], index = np.around(Echelle_R[::-1], decimals=3), columns = np.around(Echelle_s4, decimals=3))

    return data0, min, max


#===================
#To plot the figure with 9 heatmaps, for 3 values of f and 3 values of p_in, with r taking values in "Echelle_R" and s_4 taking values in "Echelle_s4", for a given value of s_3
#===================

def comparaison_heatmaps(s3, Echelle_R, Echelle_s4):

    #Divide the figure into 9 panels
    fig, axs = plt.subplots(3,3, figsize=(15, 8), sharex=True, sharey=True)
    #Position the color bar at
    cbar_ax = fig.add_axes([.95, .15, .03, .7])

    #Computation of the eigenvalues for each panel, and their min and max, used later to plot the heatmaps.

    #f=0, pin=0
    data11, min11, max11 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0, Echelle_s4, s3, f=0)

    #f=0, pin=0.5
    data12, min12, max12 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0.5, Echelle_s4, s3, f=0)

    #f=0, pin=1
    data13, min13, max13 = evol_vap_overdom_s4_r_nograph(Echelle_R, 1, Echelle_s4, s3, f=0)

    #f=0.5, pin=0
    data21, min21, max21 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0, Echelle_s4, s3, f=0.5)

    #f=0.5, pin=0.5
    data22, min22, max22 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0.5, Echelle_s4, s3, f=0.5)

    #f=0.5, pin=1
    data23, min23, max23 = evol_vap_overdom_s4_r_nograph(Echelle_R, 1, Echelle_s4, s3, f=0.5)

    #f=1, pin=0
    data31, min31, max31 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0, Echelle_s4, s3, f=1)

    #f=1, pin=0.5
    data32, min32, max32 = evol_vap_overdom_s4_r_nograph(Echelle_R, 0.5, Echelle_s4, s3, f=1)

    #f=1, pin=1
    data33, min33, max33 = evol_vap_overdom_s4_r_nograph(Echelle_R, 1, Echelle_s4, s3, f=1)

    mini = min(min11, min12, min13, min21, min22, min23, min31, min32, min33)
    maxi = max(max11, max12, max13, max21, max22, max23, max31, max32, max33)


    sns.heatmap(data11, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[0,0])
    sns.heatmap(data12, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[0,1])
    sns.heatmap(data13, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[0,2])
    sns.heatmap(data21, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[1,0])
    sns.heatmap(data22, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[1,1])
    sns.heatmap(data23, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[1,2])
    sns.heatmap(data31, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[2,0])
    sns.heatmap(data32, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[2,1])
    sns.heatmap(data33, vmin=mini, vmax=maxi, cmap='PuOr_r', center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=100, yticklabels=200, ax=axs[2,2])

    #Common xlabels and ylabels
    ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(r'Selection coefficient ($s_4$)', fontsize=20)
    plt.ylabel('Recombination rate (r)', fontsize=20)
    ax.xaxis.set_label_coords(.5, -.15)
    ax.yaxis.set_label_coords(-.08, .6)


    #Indicating the values of f and pin
    axs[2,0].set_xlabel(r'$p_{in}=0$', fontsize = 18.0)
    axs[2,1].set_xlabel(r'$p_{in}=0.5$', fontsize = 18.0)
    axs[2,2].set_xlabel(r'$p_{in}=1$', fontsize = 18.0)
    axs[0,0].set_ylabel('f=0', fontsize = 18.0)
    axs[1,0].set_ylabel('f=0.5', fontsize = 18.0)
    axs[2,0].set_ylabel('f=1', fontsize = 18.0)

    #Hide ticks for 'inside' heatmaps
    axs[0,0].tick_params(labelcolor='none', axis='x', top=False, bottom=False, left=False, right=False)
    axs[0,1].tick_params(labelcolor='none', axis='both', top=False, bottom=False, left=False, right=False)
    axs[0,2].tick_params(labelcolor='none', axis='both', top=False, bottom=False, left=False, right=False)
    axs[1,0].tick_params(labelcolor='none', axis='x', top=False, bottom=False, left=False, right=False)
    axs[1,1].tick_params(labelcolor='none', axis='both', top=False, bottom=False, left=False, right=False)
    axs[1,2].tick_params(labelcolor='none', axis='both', top=False, bottom=False, left=False, right=False)
    axs[2,1].tick_params(labelcolor='none', axis='y', top=False, bottom=False, left=False, right=False)
    axs[2,2].tick_params(labelcolor='none', axis='y', top=False, bottom=False, left=False, right=False)

    #To indicate the biological meaning of the colors alongside the colorbar
    macouleur='0.33'
    plt.arrow(1.05, 0.07, 0, 0.87, clip_on=False, color=macouleur, head_width=0.01)
    plt.text(1.02, 0.25, 'Increasing mutant maintenance', fontdict=None, clip_on=False, rotation = 90, color=macouleur, fontweight = 'semibold', fontsize = 15)
    plt.text(1.08, 1, r'$\rho$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'heavy', fontsize=15)

    #To add the line y=2x for comparison with previous references
    #To add a line over a heatmap, need to specify the range of this line in terms of position on the axes, not values of datas. So to have the line y=2x, since x-axis is [0;0.1] and y-axis [0;0.5], we need to plot a line that covers [0;2/5] of the graph. Here, the origin is at the top left (I don't know why), thus the coordinates need to be (1, 3/5) rather than (0,2/5). The same system of coordinates is used to place the text.
    n = len(Echelle_s4)
    m = len(Echelle_R)
    Y = np.linspace(n+1, m*3/5, n+1)
    #To add the "r=2s4" on top of the line, at the right position.
    echelle_x = n/100
    echelle_y = m/100

    #Plotting the line y=2x on the six bottom panels
    for k in range(3):
        axs[2,k].plot(Y, color=macouleur, linewidth=3, linestyle='dashed')
        axs[2,k].text(70*echelle_x, 67*echelle_y, r'$r=2s_4$', fontdict=None, clip_on=False, rotation = 14, color=macouleur, fontweight = 'semibold', fontsize = 15)
        axs[1,k].plot(Y, color=macouleur, linewidth=3, linestyle='dashed')
        axs[1,k].text(70*echelle_x, 67*echelle_y, r'$r=2s_4$', fontdict=None, clip_on=False, rotation = 14, color=macouleur, fontweight = 'semibold', fontsize = 15)


    plt.savefig("Fig2_Heatmaps_VapOverdom.png", bbox_inches="tight")
    plt.show()

#===================
#Parameters and commands
#===================

s3 = 0.1
Ech_R = np.linspace(10**(-10), 0.5, 1000)
Ech_s4 = np.linspace(0, s3, 1000)

comparaison_heatmaps(s3, Ech_R, Ech_s4)
