import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd


#File to generate figures 4 and S3 : bar charts of the time to extinction of the branching process in the partial dominance case.
#The overdominance case is treated in a separate file as the legends and number of simultaneous bar charts are not always the same. However, the functions used are the same for the most part.

#To plot histograms for 3 values of r on the same graph
def histolog_multiples_r(Tableau1, Tableau2, Tableau3):

    #Extracting parameters values
    r1, pin, s, h, f, X01, X02, nb_traj = Tableau1[0:8]
    r2 = Tableau2[0]
    r3 = Tableau3[0]

    #Extracting the times to build a histogram with
    Temps1 = Tableau1[8:]
    Temps2 = Tableau2[8:]
    Temps3 = Tableau3[8:]


    #To plot indicative vertical lines
    maxi1 = np.max(Temps1)
    maxi2 = np.max(Temps2)
    maxi3 = np.max(Temps3)
    print(maxi1, maxi2, maxi3)
    moy1 = np.mean(Temps1)
    moy2 = np.mean(Temps2)
    moy3 = np.mean(Temps3)
    print(moy1, moy2, moy3)
    q75_1 = np.percentile(Temps1, 75)
    q75_2 = np.percentile(Temps2, 75)
    q75_3 = np.percentile(Temps3, 75)
    q99_1 = np.percentile(Temps1, 99)
    q99_2 = np.percentile(Temps2, 99)
    q99_3 = np.percentile(Temps3, 99)

    #As extinction times span several order of magnitude, x-axis is on a log scale
    logintervalles = np.logspace(0, 6, 100, base = 10)

    #Plotting the bar chart
    couleurs = ['deepskyblue', 'coral', 'limegreen']
    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(111)
    ax.hist(Temps1, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.3f$' %r1, color = couleurs[0])
    ax.hist(Temps2, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.1f$' %r2, color = couleurs[1])
    ax.hist(Temps3, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.1f$' %r3, color=couleurs[2])
    ax.set_xscale('log')

    ax.set_xlabel('Purging time of the deleterious allele', fontsize=15)
    ax.set_ylabel('Density', fontsize=15)

    #Plotting the lines
    ax.vlines(q75_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q75 ($r=%.3f$)' %r1, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    ax.vlines(q75_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q75 ($r=%.1f$)' %r2, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    ax.vlines(q75_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q75 ($r=%.1f$)' %r3, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    ax.vlines(moy1, ymin=0, ymax=0.175, color= couleurs[0], label = r'mean ($r=%.3f$)' %r1, linewidth = 5)
    ax.vlines(moy2, ymin=0, ymax=0.175, color= couleurs[1], label = r'mean ($r=%.1f$)' %r2, linewidth = 5)
    ax.vlines(moy3, ymin=0, ymax=0.175, color= couleurs[2], label = r'mean ($r=%.1f$)' %r3, linewidth = 5)
    ax.vlines(q99_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q99 ($r=%.3f$)' %r1, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    ax.vlines(q99_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q99 ($r=%.1f$)' %r2, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    ax.vlines(q99_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q99 ($r=%.1f$)' %r3, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    ax.vlines(maxi1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'max ($r=%.3f$)' %r1, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    ax.vlines(maxi2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'max ($r=%.1f$)' %r2, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    ax.vlines(maxi3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'max ($r=%.1f$)' %r3, linestyles = 'dotted', linewidth = 5, alpha=0.7)


    ax.legend(prop={'size':15},bbox_to_anchor=(0.5, 1.24), loc='upper center', ncol=5)
    fig.subplots_adjust(top=0.8)
    fig.show()
    fig.savefig("Fig4_Histo_contresel_r_s01.png")



#To plot histograms for 3 values of f on the same graph
def histolog_multiples_f(Tableau1, Tableau2, Tableau3):

    #Extracting parameters values
    r, pin, s, h, f1, X01, X02, nb_traj = Tableau1[0:8]
    f2 = Tableau2[4]
    f3 = Tableau3[4]

    #Extracting the times to build a histogram with
    Temps1 = Tableau1[8:]
    Temps2 = Tableau2[8:]
    Temps3 = Tableau3[8:]


    #To plot indicative vertical lines
    maxi1 = np.max(Temps1)
    maxi2 = np.max(Temps2)
    maxi3 = np.max(Temps3)
    print(maxi1, maxi2, maxi3)
    moy1 = np.mean(Temps1)
    moy2 = np.mean(Temps2)
    moy3 = np.mean(Temps3)
    print(moy1, moy2, moy3)
    q75_1 = np.percentile(Temps1, 75)
    q75_2 = np.percentile(Temps2, 75)
    q75_3 = np.percentile(Temps3, 75)
    q99_1 = np.percentile(Temps1, 99)
    q99_2 = np.percentile(Temps2, 99)
    q99_3 = np.percentile(Temps3, 99)

    #As extinction times span several order of magnitude, x-axis is on a log scale
    logintervalles = np.logspace(0, 6, 100, base = 10)

    #Plotting the bar chart
    couleurs = [ 'deepskyblue', 'coral', 'limegreen']
    plt.figure(figsize=(18,10))
    plt.hist(Temps1, bins=logintervalles, density=1, alpha = 1, label=r'$f=%.d$' %f1, color = couleurs[0])
    plt.hist(Temps2, bins=logintervalles, density=1, alpha = 1, label=r'$f=%.1f$' %f2, color = couleurs[1])
    plt.hist(Temps3, bins=logintervalles, density=1, alpha = 1, label=r'$f=%.d$' %f3, color=couleurs[2])
    plt.xscale('log')

    plt.xlabel('Purging time of the deleterious allele', fontsize=15)
    plt.ylabel('Density', fontsize=15)

    #Plotting the vertical lines
    plt.vlines(q75_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q75 ($f=%.3f$)' %f1, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q75 ($f=%.1f$)' %f2, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q75 ($f=%.1f$)' %f3, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(moy1, ymin=0, ymax=0.175, color= couleurs[0], label = r'mean ($f=%.3f$)' %f1, linewidth = 5)
    plt.vlines(moy2, ymin=0, ymax=0.175, color= couleurs[1], label = r'mean ($f=%.1f$)' %f2, linewidth = 5)
    plt.vlines(moy3, ymin=0, ymax=0.175, color= couleurs[2], label = r'mean ($f=%.1f$)' %f3, linewidth = 5)
    plt.vlines(q99_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q99 ($f=%.3f$)' %f1, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q99 ($f=%.1f$)' %f2, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q99 ($f=%.1f$)' %f3, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(maxi1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'max ($f=%.3f$)' %f1, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'max ($f=%.1f$)' %f2, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'max ($f=%.1f$)' %f3, linestyles = 'dotted', linewidth = 5, alpha=0.7)


    plt.legend(prop={'size':15},bbox_to_anchor=(0.5, 1.24), loc='upper center', ncol=5)
    plt.subplots_adjust(top=0.8)
    plt.show()
    plt.savefig("FigS3_Histo_contresel_f_s01.png")


#==================================
#Commands
#==================================
#Commands used to display figures. Here the imported tables are the ones I used to plot the figures. Other set of data can be obtained by launching the SimuBranchingProcess_PartialDominance file, with wished parameters values. Each launch of the latter file will give a different set of data, as the process is stochastic.

#Figure 4
T1 = np.loadtxt('Thisto_r0001_s01', delimiter = ';')
T2 = np.loadtxt('Thisto_r01_s01', delimiter = ';')
T3 = np.loadtxt('Thisto_r05_s01', delimiter = ';')
#histolog_multiples_r(T1, T2, T3)

#Figure S3
T4 = np.loadtxt('Thisto_f0_s01', delimiter = ';')
T5 = np.loadtxt('Thisto_f1_s01', delimiter = ';')
T6 = np.loadtxt('Thisto_f05_s01', delimiter = ';')
histolog_multiples_f(T5, T6, T4)
