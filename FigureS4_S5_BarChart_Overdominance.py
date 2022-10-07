import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd

#File to generate figures S4 and S5 : bar charts of the time to extinction of the branching process in the overdominance case.
#The partial dominance case is treated in a separate file as the legends and number of simultaneous bar charts are not always the same. However, the functions used are the same for the most part.

#To plot histograms for 3 sets of (r, s4) values on the same graph
def histolog_multiples_r(Tableau1, Tableau2, Tableau3):

    #Extracting parameters values
    r1, pin, s3, s4, f, X01, X02, nb_traj = Tableau1[0:8]
    r2, s4_2 = Tableau2[0], Tableau2[3]
    r3, s4_3 = Tableau3[0], Tableau3[3]

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
    plt.hist(Temps1, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.1f$, $s_4 = %.3f$' %(r1, s4), color = couleurs[0])
    plt.hist(Temps2, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.1f$, $s_4 = %.2f$' %(r2, s4_2), color = couleurs[1])
    plt.hist(Temps3, bins=logintervalles, density=1, alpha = 1, label=r'$r=%.1f$, $s_4 = %.3f$' %(r3, s4_3), color=couleurs[2])
    plt.xscale('log')

    plt.xlabel('Purging time of the deleterious allele', fontsize=15)
    plt.ylabel('Density', fontsize=15)

    #Plotting the vertical lines
    plt.vlines(q75_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q75 ($r=%.1f$, $s_4=%.3f$)' %(r1, s4), linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q75 ($r=%.1f$, $s_4=%.2f$)' %(r2, s4_2), linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q75 ($r=%.1f$, $s_4=%.3f$)' %(r3, s4_3), linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(moy1, ymin=0, ymax=0.175, color= couleurs[0], label = r'mean ($r=%.1f$, $s_4=%.3f$)' %(r1, s4), linewidth = 5)
    plt.vlines(moy2, ymin=0, ymax=0.175, color= couleurs[1], label = r'mean ($r=%.1f$, $s_4=%.2f$)' %(r2, s4_2), linewidth = 5)
    plt.vlines(moy3, ymin=0, ymax=0.175, color= couleurs[2], label = r'mean ($r=%.1f$, $s_4=%.3f$)' %(r3, s4_3), linewidth = 5)
    plt.vlines(q99_1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'q99 ($r=%.1f$, $s_4=%.3f$)' %(r1, s4), linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'q99 ($r=%.1f$, $s_4=%.2f$)' %(r2, s4_2), linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'q99 ($r=%.1f$, $s_4=%.3f$)' %(r3, s4_3), linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(maxi1, ymin = 0, ymax = 0.175, color = couleurs[0], label = r'max ($r=%.1f$, $s_4=%.3f$)' %(r1, s4), linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi2, ymin = 0, ymax = 0.175, color = couleurs[1], label = r'max ($r=%.1f$, $s_4=%.2f$)' %(r2, s4_2), linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi3, ymin = 0, ymax = 0.175, color = couleurs[2], label = r'max ($r=%.1f$, $s_4=%.3f$)' %(r3, s4_3), linestyles = 'dotted', linewidth = 5, alpha=0.7)


    plt.legend(prop={'size':15},bbox_to_anchor=(0.5, 1.24), loc='upper center', ncol=5)
    plt.subplots_adjust(top=0.8)
    plt.show()
    plt.savefig("FigS4_Histo_overdom_r.png")



#To plot histograms for 4 values of f on the same graph
def histolog_multiples_f(Tableau1, Tableau2, Tableau3, Tableau4):

    #Extracting parameters values
    r, pin, s3, s4, f1, X01, X02, nb_traj = Tableau1[0:8]
    f2 = Tableau2[4]
    f3 = Tableau3[4]
    f4 = Tableau4[4]

    #Extracting the times to build a histogram with
    Temps1 = Tableau1[8:]
    Temps2 = Tableau2[8:]
    Temps3 = Tableau3[8:]
    Temps4 = Tableau4[8:]

    #To plot indicative vertical lines
    maxi1 = np.max(Temps1)
    maxi2 = np.max(Temps2)
    maxi3 = np.max(Temps3)
    maxi4 = np.max(Temps4)
    print(maxi1, maxi2, maxi3, maxi4)
    moy1 = np.mean(Temps1)
    moy2 = np.mean(Temps2)
    moy3 = np.mean(Temps3)
    moy4 = np.mean(Temps4)
    print(moy1, moy2, moy3, moy4)
    q75_1 = np.percentile(Temps1, 75)
    q75_2 = np.percentile(Temps2, 75)
    q75_3 = np.percentile(Temps3, 75)
    q75_4 = np.percentile(Temps4, 75)
    q99_1 = np.percentile(Temps1, 99)
    q99_2 = np.percentile(Temps2, 99)
    q99_3 = np.percentile(Temps3, 99)
    q99_4 = np.percentile(Temps4, 99)


    #As extinction times span several order of magnitude, x-axis is on a log scale
    logintervalles = np.logspace(0, 6, 100, base = 10)

    #Plotting the bar chart
    couleurs = ['limegreen', 'deepskyblue', 'coral', 'gold']
    plt.figure(figsize=(18,10))
    plt.hist(Temps1, bins=logintervalles, density=1, alpha = 0.5, label=r'$f=%.d$' %f1, color = couleurs[0])
    plt.hist(Temps2, bins=logintervalles, density=1, alpha = 0.7, label=r'$f=%.1f$' %f2, color = couleurs[1])
    plt.hist(Temps3, bins=logintervalles, density=1, alpha = 0.9, label=r'$f=%.1f$' %f3, color=couleurs[2])
    plt.hist(Temps4, bins=logintervalles, density=1, alpha = 0.7, label=r'$f=%.1f$' %f4, color=couleurs[3])
    plt.xscale('log')


    plt.xlabel('Purging time of the deleterious allele', fontsize=15)
    plt.ylabel('Density', fontsize=15)


    #Plotting the vertical lines
    plt.vlines(q75_1, ymin = 0, ymax = 0.2, color = couleurs[0], label = r'q75 ($f=%.d$)' %f1, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_2, ymin = 0, ymax = 0.2, color = couleurs[1], label = r'q75 ($f=%.1f$)' %f2, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_3, ymin = 0, ymax = 0.2, color = couleurs[2], label = r'q75 ($f=%.1f$)' %f3, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(q75_4, ymin = 0, ymax = 0.2, color = couleurs[3], label = r'q75 ($f=%.1f$)' %f4, linestyles = (0, (2, 4)), linewidth = 4, alpha=0.7)
    plt.vlines(moy1, ymin=0, ymax=0.2, color= couleurs[0], label = r'mean ($f=%.d$)' %f1, linewidth = 5)
    plt.vlines(moy2, ymin=0, ymax=0.2, color= couleurs[1], label = r'mean ($f=%.1f$)' %f2, linewidth = 5)
    plt.vlines(moy3, ymin=0, ymax=0.2, color= couleurs[2], label = r'mean ($f=%.1f$)' %f3, linewidth = 5)
    plt.vlines(moy4, ymin=0, ymax=0.2, color= couleurs[3], label = r'mean ($f=%.1f$)' %f4, linewidth = 5)
    plt.vlines(q99_1, ymin = 0, ymax = 0.2, color = couleurs[0], label = r'q99 ($f=%.3f$)' %f1, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_2, ymin = 0, ymax = 0.2, color = couleurs[1], label = r'q99 ($f=%.1f$)' %f2, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_3, ymin = 0, ymax = 0.2, color = couleurs[2], label = r'q99 ($f=%.1f$)' %f3, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(q99_4, ymin = 0, ymax = 0.2, color = couleurs[2], label = r'q99 ($f=%.1f$)' %f4, linestyles = 'dashed', linewidth = 5, alpha=0.7)
    plt.vlines(maxi1, ymin = 0, ymax = 0.2, color = couleurs[0], label = r'max ($f=%.d$)' %f1, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi2, ymin = 0, ymax = 0.2, color = couleurs[1], label = r'max ($f=%.1f$)' %f2, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi3, ymin = 0, ymax = 0.2, color = couleurs[2], label = r'max ($f=%.1f$)' %f3, linestyles = 'dotted', linewidth = 5, alpha=0.7)
    plt.vlines(maxi4, ymin = 0, ymax = 0.2, color = couleurs[3], label = r'max ($f=%.1f$)' %f4, linestyles = 'dotted', linewidth = 5, alpha=0.7)


    plt.legend(prop={'size':15},bbox_to_anchor=(0.5, 1.26), loc='upper center', ncol=5)
    plt.subplots_adjust(top=0.8)
    plt.show()
    plt.savefig("FigS5_Histo_overdom_f.png")

#==================================
#Commands
#==================================


#Figure S4
T1 = np.loadtxt('T_histo_overdom_r01_s40001', delimiter=';')
T2 = np.loadtxt('T_histo_overdom_r05_s4001', delimiter=';')
T3 = np.loadtxt('T_histo_overdom_r05_s40001', delimiter=';')
#histolog_multiples_r(T1, T2, T3)

#Figure S5
T4 = np.loadtxt('Thisto_overdom_f01', delimiter=';')
T5= np.loadtxt('Thisto_overdom_f05', delimiter=';')
T6 = np.loadtxt('Thisto_overdom_f09', delimiter=';')
T7 = np.loadtxt('Thisto_overdom_f1', delimiter=';')
#histolog_multiples_f(T7, T6, T5, T4)
