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
#Figure 3 : Plot the heatmaps on the same scale for the two scenarios, for varing f and pin
#===============================================================================================

def comparaison_heatmap_derivee(Echelle_pin, Echelle_f, r, h, s, s4, s3):

    #Computing the derivative values for both heatmaps
    min_contresel, max_contresel, les_deriv_contresel = evol_deriv_rho_contresel(Echelle_pin, Echelle_f, r, h, s)
    min_overdom, max_overdom, les_deriv_overdom = evol_deriv_rho_overdom(Echelle_pin, Echelle_f, r, s4, s3)

    #Rescaling the data to have all values in [0,1] in order to compare the two heatmaps on the same scale
    deriv_contresel_rescale = les_deriv_contresel/min_contresel
    deriv_overdom_rescale = les_deriv_overdom/min_overdom

    #Editing datas so that it can be understood by seaborn
    data_contresel = pd.DataFrame(deriv_contresel_rescale, columns = np.around(Echelle_f, decimals=3), index = np.around(Echelle_pin[::-1], decimals=3))
    data_overdom = pd.DataFrame(deriv_overdom_rescale, index = np.around(Echelle_pin[::-1], decimals=3), columns = np.around(Echelle_f, decimals=3))

    #Computing the new min and max to set bounds on the colorbar
    mini = min(max_contresel/min_contresel, max_overdom/min_overdom)
    maxi = 1

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
    plt.text(0.47, 1, r'$\frac{\partial \rho}{\partial r}|_{r=0.5}$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'heavy', fontsize=20)

    plt.show()
    plt.savefig("Fig3_derivee_contresel_overdom.png", bbox_inches="tight")


#===============================================================================================
#Figure S1 : Plot the heatmaps of the dominant eigenvalue derivative for the partial dominance scenario, for varying f and h, and varying f and s
#===============================================================================================


#Computing the derivative for the heatmap f/h
def evol_deriv_rho_contresel_h(pin, Echelle_f, r, Echelle_h, s):

    n = len(Echelle_f)
    m = len(Echelle_h)
    les_deriv = np.zeros((m,n))

    for k in range(m):
        h = Echelle_h[m-1-k]#So that the y-axis is in the right order (increasing h when going up)
        for l in range(n):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_contresel(r, pin, h, s, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv

#Computing the derivative for the heatmap f/s
def evol_deriv_rho_contresel_s(pin, Echelle_f, r, h, Echelle_s):

    n = len(Echelle_f)
    m = len(Echelle_s)
    les_deriv = np.zeros((m,n))

    for k in range(m):
        s = Echelle_s[m-1-k]#So that the y-axis is in the right order (increasing s when going up)
        for l in range(n):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_contresel(r, pin, h, s, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv

#Plotting the heatmaps for figure S1 : left panel f/h and right panel f/s
def comparaison_heatmap_derivee_contresel_supp():

    Ech_h = np.linspace(0, 0.5, 1000)
    Ech_f = np.linspace(0, 1, 1001)
    Ech_s = np.linspace(0, 0.1, 1000)
    r = 0.5
    h = 0.3
    pin = 0.5
    s = 0.1

    #Computing the derivative values for both heatmaps
    min_h, max_h, les_deriv_h = evol_deriv_rho_contresel_h(pin, Ech_f, r, Ech_h, s)
    min_s, max_s, les_deriv_s = evol_deriv_rho_contresel_s(pin, Ech_f, r, h, Ech_s)

    #Rescaling the data to have all values in [0,1] in order to compare the two heatmaps on the same scale
    deriv_h_rescale = les_deriv_h/min_h
    deriv_s_rescale = les_deriv_s/min_s

    #Editing datas so that it can be understood by seaborn
    data_h = pd.DataFrame(deriv_h_rescale, columns = np.around(Ech_f, decimals=3), index = np.around(Ech_h[::-1], decimals=3))
    data_s = pd.DataFrame(deriv_s_rescale, index = np.around(Ech_s[::-1], decimals=3), columns = np.around(Ech_f, decimals=3))

    #Computing the new min and max to set bounds on the colorbar
    mini = min(max_h/min_h, max_s/min_s)
    maxi = 1

    #Editing figure to have two panels and a single colorbar
    fig, axs = plt.subplots(1,2, figsize=(18, 10), sharex=True, sharey=False)
    fig.subplots_adjust(wspace=0.6)
    cbar_ax = fig.add_axes([.48, .15, .03, .7])

    #Plotting the first heatmap (f/h)
    sns.heatmap(data_h, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[0])
    axs[0].set_title(r'$s = %.1f, p_{in} = %.1f$' %(s, pin), fontsize=20)
    axs[0].set_ylabel('Dominance coefficient (h)', fontsize=20)

    #Plotting the second heatmap (f/s)
    sns.heatmap(data_s, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[1])
    axs[1].set_title(r'$p_{in} = %.1f, h =%.1f$' %(pin, h), fontsize = 20)
    axs[1].set_ylabel(r'Selection coefficient ($s$)', fontsize=20)

    #Common xlabels and ylabels
    ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Selfing rate (f)', fontsize=20)
    ax.xaxis.set_label_coords(.5, -.11)
    ax.yaxis.set_label_coords(-.08, .5)

    #To indicate the biological meaning of the colors alongside the colorbar
    plt.arrow(0.44, 0.07, 0, 0.87, clip_on=False, color='dimgrey', head_width=0.01)
    plt.text(0.42, 0.25, 'Increasing mutant maintenance', fontdict=None, clip_on=False, rotation = 90, color='dimgrey', fontweight = 'semibold', fontsize = 15)
    plt.text(0.44, 1, r'$\frac{\partial \rho}{\partial r}|_{r=0.5}$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'heavy', fontsize=20)

    plt.show()
    plt.savefig("FigS1_derivee_contresel_h_s.png", bbox_inches="tight")



#===============================================================================================
#Figure S2 : Plot the heatmaps of the dominant eigenvalue derivative for the overdominance scenario, for varying f and s3, and varying f and s4
#===============================================================================================

#Computing the derivative for the heatmap f/s3
def evol_deriv_rho_overdom_s3(pin, Echelle_f, r, s4, Echelle_s3):

    m = len(Echelle_f)
    n = len(Echelle_s3)
    les_deriv = np.zeros((n,m))

    for k in range(n):
        s3 = Echelle_s3[n-1-k]#Pour que l'axe des ordonnées soit dans le sens croissant en "montant"
        for l in range(m):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_overdom(r, pin, s3, s4, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv

#Computing the derivative for the heatmap f/s4
def evol_deriv_rho_overdom_s4(pin, Echelle_f, r, Echelle_s4, s3):

    m = len(Echelle_f)
    n = len(Echelle_s4)
    les_deriv = np.zeros((n,m))

    for k in range(n):
        s4 = Echelle_s4[n-1-k]#Pour que l'axe des ordonnées soit dans le sens croissant en "montant"
        for l in range(m):
            f = Echelle_f[l]
            calcul_deriv = derivee_rho_overdom(r, pin, s3, s4, f)
            les_deriv[k][l] = calcul_deriv

    min = np.min(les_deriv)
    max = np.max(les_deriv)


    return min, max, les_deriv

#Plotting the heatmaps for figure S2 : left panel f/s3 and right panel f/s4
def comparaison_heatmap_derivee_overdom_supp():

    #Careful to keep s3 > s4
    Ech_s3 = np.linspace(0, 0.1, 1000)
    Ech_f = np.linspace(0, 1, 1001)
    Ech_s4 = np.linspace(0, 0.1, 1000)
    r = 0.5
    s3 = 0.1
    pin = 0.5
    s4 = 0.001

    #Computing the derivative values for both heatmaps
    min_s3, max_s3, les_deriv_s3 = evol_deriv_rho_overdom_s3(pin, Ech_f, r, s4, Ech_s3)
    min_s4, max_s4, les_deriv_s4 = evol_deriv_rho_overdom_s4(pin, Ech_f, r, Ech_s4, s3)

    #Rescaling the data to have all values in [0,1] in order to compare the two heatmaps on the same scale
    deriv_s3_rescale = les_deriv_s3/min_s3
    deriv_s4_rescale = les_deriv_s4/min_s4

    #Editing datas so that it can be understood by seaborn
    data_s3 = pd.DataFrame(deriv_s3_rescale, columns = np.around(Ech_f, decimals=3), index = np.around(Ech_s3[::-1], decimals=3))
    data_s4 = pd.DataFrame(deriv_s4_rescale, index = np.around(Ech_s4[::-1], decimals=3), columns = np.around(Ech_f, decimals=3))

    #Computing the new min and max to set bounds on the colorbar
    mini = min(max_s3/min_s3, max_s4/min_s4)
    maxi = 1

    #Editing figure to have two panels and a single colorbar
    fig, axs = plt.subplots(1,2, figsize=(18, 10), sharex=False, sharey=False)
    fig.subplots_adjust(wspace=0.6)
    cbar_ax = fig.add_axes([.48, .15, .03, .7])

    #Plotting the first heatmap (f/s3)
    sns.heatmap(data_s3, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[0])
    axs[0].set_title(r'$s_4 = %.3f$, $p_{in} = %.1f$' %(s4, pin), fontsize=20)
    axs[0].set_ylabel(r'$bb$ selection coefficient ($s_3$)', fontsize=20)

    #Plotting the first heatmap (f/s4)
    sns.heatmap(data_s4, vmin=mini, vmax=maxi, cmap="PuOr_r", center = 0, cbar=True, cbar_kws=None, cbar_ax=cbar_ax, square=False, xticklabels=200, yticklabels=200, ax=axs[1])
    axs[1].set_title(r'$p_{in} = %.1f$, $s_3 =%.1f$' %(pin, s3), fontsize=20)
    axs[1].set_ylabel(r'$BB$ selection coefficient ($s_4$)', fontsize=20)


    #Common xlabels and ylabels
    ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Selfing rate (f)', fontsize=20)
    ax.xaxis.set_label_coords(.5, -.11)
    ax.yaxis.set_label_coords(-.08, .5)

    #To indicate the biological meaning of the colors alongside the colorbar
    plt.arrow(0.44, 0.07, 0, 0.87, clip_on=False, color='dimgrey', head_width=0.01)
    plt.text(0.42, 0.25, 'Increasing mutant maintenance', fontdict=None, clip_on=False, rotation = 90, color='dimgrey', fontweight = 'semibold', fontsize = 15)
    plt.text(0.44, 1, r'$\frac{\partial \rho}{\partial r}|_{r=0.5}$', fontdict=None, clip_on=False, rotation = 0, color='black', fontweight = 'heavy', fontsize=20)



    plt.show()
    plt.savefig("FigS2_derivee_overdom_s3_s4.png", bbox_inches="tight")



#===============================================================================================
#Parameters and Commands
#===============================================================================================
Ech_pin = np.linspace(0,1, 1001)
Ech_f = np.linspace(0, 1, 1001)
r = 0.5
h = 0.3
s = 0.1
s4 = 0.05
s3 = 0.1

#Remark : for figures S1 and S2, parameters are included in the function, whereas they have to be specified for figure 3

#Figure 3
#comparaison_heatmap_derivee(Ech_pin, Ech_f, r, h, s, s4, s3)

#Figure S1
#comparaison_heatmap_derivee_contresel_supp()

#Figure S2
#comparaison_heatmap_derivee_overdom_supp()




