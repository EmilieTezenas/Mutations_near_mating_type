import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import math
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

#In this file, "overdom" stands for the overdominance scenario, and "contresel" for the partial dominance scenario

#This file computes and shows the difference between the derivative at r=0 and the derivative at r=0.5

#===============================================================================================
#Compute the derivative of the dominant eigenvalue for the Permanent Disadvantage scenario, and its variation along a pin-scale and a f-scale to plot the heatmap
#===============================================================================================

def derivee_rho_contresel(r, pin, h, s, f):

    pout = 1-pin
    alpha = 2 - r - pin*(1-r)

    beta = f*(-r*(1-h*s)*alpha + 2*(1-s)) - 2*(1+h*s)
    beta_deriv = -f*(1-h*s)*(1+(1-pin)*(1-2*r))

    delta = (beta + 4*h*s)**2 - 8*f*s*r*alpha*(1-h)*(1-h*s)
    delta_deriv = 2*beta_deriv*(beta+4*h*s) - 8*f*s*(1-h)*(1-h*s)*(1+(1-pin)*(1-2*r))

    if delta == 0:
        lambda_deriv = 0
    else :
        lambda_deriv = 0.25*(beta_deriv + delta_deriv/(2*np.sqrt(delta)))

    return lambda_deriv

#Gives the data for the heatmap + min and max

def evol_deriv_rho_contresel(Echelle_pin, Echelle_f, r, h, s):

    n = len(Echelle_f)
    m = len(Echelle_pin)
    les_deriv = np.zeros((m,n))

    for k in range(m):
        pin = Echelle_pin[m-1-k]#Pour que l'axe des ordonnées soit dans le sens croissant en "montant"
        for l in range(n):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_contresel(r, pin, h, s, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv


#===============================================================================================
#Compute the derivative of the dominant eigenvalue for the Permanent Disadvantage scenario, and its variation along a pin-scale and a f-scale to plot the heatmap
#===============================================================================================

def derivee_rho_overdom(r, pin, s3, s4, f):

    alpha= pin*(1-r)+r-2
    beta = f*(r*alpha - 2*s3 +2) - 2 + 4*s4
    deriv_beta = 2*f*(1-pin)*r + f*(pin-2)

    Delta = f**2*(alpha*r - 2*s3 +2)**2 + 4*f*(2*s3 - 2 + r*(2*s3-1)*alpha) + 4
    deriv_Delta = 2*deriv_beta*(beta-4*s4) + 8*f*s3*(2*(1-pin)*r + pin-2)

    if Delta == 0:
        deriv_lambda = 0
    else :
        deriv_lambda = 0.25*(deriv_beta + 0.5*deriv_Delta/np.sqrt(Delta))

    return deriv_lambda

#Gives the data for the heatmap + min and max

def evol_deriv_rho_overdom(Echelle_pin, Echelle_f, r, s4, s3):

    m = len(Echelle_f)
    n = len(Echelle_pin)
    les_deriv = np.zeros((n,m))

    for k in range(n):
        pin = Echelle_pin[n-1-k]#Pour que l'axe des ordonnées soit dans le sens croissant en "montant"
        for l in range(m):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_overdom(r, pin, s3, s4, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv

#===============================================================================================
#Figure S3 : Plot the difference between the derivatives at r=0.5 and r=0, for the two scenarios, for varing f and pin
#===============================================================================================

def comparaison_heatmap_derivee(Echelle_pin, Echelle_f, r, h, s, s4, s3):


    #Computing the derivative at r=0 values for both heatmaps
    min0_contresel, max0_contresel, les_deriv0_contresel = evol_deriv_rho_contresel(Echelle_pin, Echelle_f, 0, h, s)
    min0_overdom, max0_overdom, les_deriv0_overdom = evol_deriv_rho_overdom(Echelle_pin, Echelle_f, 0, s4, s3)

    #Computing the derivative at r=0.5 values for both heatmaps
    min05_contresel, max05_contresel, les_deriv05_contresel = evol_deriv_rho_contresel(Echelle_pin, Echelle_f, 0.5, h, s)
    min05_overdom, max05_overdom, les_deriv05_overdom = evol_deriv_rho_overdom(Echelle_pin, Echelle_f, 0.5, s4, s3)

    #Computing the difference
    deriv_contresel_diff = les_deriv05_contresel - les_deriv0_contresel
    deriv_overdom_diff = les_deriv05_overdom - les_deriv0_overdom

    #Editing datas so that it can be understood by seaborn
    data_contresel = pd.DataFrame(deriv_contresel_diff, columns = np.around(Echelle_f, decimals=3), index = np.around(Echelle_pin[::-1], decimals=3))
    data_overdom = pd.DataFrame(deriv_overdom_diff, index = np.around(Echelle_pin[::-1], decimals=3), columns = np.around(Echelle_f, decimals=3))

    #Computing the new min and max to set bounds on the colorbar
    mini = min(np.min(deriv_contresel_diff), np.min(deriv_overdom_diff))
    maxi = max(np.max(deriv_contresel_diff), np.max(deriv_overdom_diff))


    #Editing figure to have two panels and a single colorbar
    fig, axs = plt.subplots(1,2, figsize=(18, 10), sharex=True, sharey=False)
    fig.subplots_adjust(wspace=0.4)
    cbar_ax = fig.add_axes([.5, .15, .03, .7])

    #Plotting the first heatmap (partial dominance)
    sns.heatmap(data_contresel, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[0])
    axs[0].set_title('Partial Dominance. s = %.1f, h = %.1f' %(s, h), fontsize = 20)
    axs[0].set_yticks([0, 201, 401, 601, 801])

    #Plotting the second heatmap (overdominance)
    sns.heatmap(data_overdom, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[1])
    axs[1].set_title(r'Overdominance. $s_3$ = %.1f, $s_4$ =%.2f' %(s3, s4), fontsize = 20)
    axs[1].set_yticks([0, 201, 401, 601, 801])


    #Common xlabels and ylabels
    ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Selfing rate (f)', fontsize=20)
    plt.ylabel(r'Intratetrad selfing coefficient ($p_{in}$)', fontsize=20)
    ax.xaxis.set_label_coords(.5, -.11)
    ax.yaxis.set_label_coords(-.08, .5)

    #To indicate the biological meaning of the colors alongside the colorbar
    plt.arrow(0.47, 0.07, 0, 0.87, clip_on=False, color='dimgrey', head_width=0.01)
    plt.text(0.45, 0.25, 'Increasing mutant maintenance', fontdict=None, clip_on=False, rotation = 90, color='dimgrey', fontweight = 'semibold', fontsize = 15)
    plt.text(0.48, 1, r'$\Delta\left(\frac{\partial \rho}{\partial r}\right)$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'heavy', fontsize=20)

    plt.show()
    plt.savefig("FigS3_derivee_contresel_overdom_diff.png", bbox_inches="tight")




#===============================================================================================
#Parameters and Commands
#===============================================================================================
Ech_pin = np.linspace(0,1, 1001)
Ech_f = np.linspace(0, 1, 1001)
r = 0
h = 0.3
s = 0.1
s4 = 0.05
s3 = 0.1

#Remark : for figures S1 and S2, parameters are included in the function, whereas they have to be specified for figure 3

#Figure S3
#comparaison_heatmap_derivee(Ech_pin, Ech_f, r, h, s, s4, s3)




