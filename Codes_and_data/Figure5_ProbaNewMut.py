import numpy as np
import matplotlib.pyplot as plt

#To plot the probability that a new mutations arises before the first one is purged, for a given length of the region on which the new mutation can appear, and various recombination rates between the first mutation and the mating-type locus.


plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)


#Compute the proportion of long-enough runs, for a given set of runs
def proportion_nouvelle_mut(taux_mutation_bp, distance, distri_sauts):

    #Rate of mutation per reproduction event
    mu_ev = taux_mutation_bp*distance

    #Mean number of reproduction event to wait before seeing another mutation appear.
    nb_ev_newmut = 1/mu_ev

    longues_traj = distri_sauts[distri_sauts >= nb_ev_newmut]

    #And then return the numer of long-enough runs divided by the total number of runs.
    return len(longues_traj)/len(distri_sauts)


#Importing the table containing the runs for each set of parameter.
Tglobal = np.loadtxt('Tableau_concatene', delimiter=';')

#Organization of Tglobal :
#- Blocks of p_in (p_in = 0.5 puis p_in=0.9)
#- Then blocks of (0, 0.5, 0.8, 1)
#- Then each line for a recombination rate r (0.001, 0.01, 0.1, 0.5)
# Each line corresponds to a set of parameter, for which 100 000 independant runs were run.

#Compute the probability of a new mutation before the first is purged for each line of a given table
def calcul_proba_nouvelle_mut(taux_mut_bp, Tableau_distri, distance):

    n = len(Tableau_distri)
    les_probas=np.zeros(n)
    F = Tableau_distri[:, 4]
    Pin = Tableau_distri[:,1]
    s = Tableau_distri[0,2]
    h = Tableau_distri[0,3]
    nb_traj = Tableau_distri[0,7]


    for i in range(n):
        les_probas[i]=proportion_nouvelle_mut(taux_mut_bp, distance, Tableau_distri[i,8:])



    return les_probas, F, Pin


#Plots the probability that a new mutation arises before the first one is purged, for various f and pin, and along the recombination rate.
#ATTENTION, the program needs to be changed if the number of r, pin and f changes (if another table than Tglobal is used)
#For Tglobal : nb_R = 4, nb_pin = 2, nb_f=4
def trace_probanewmut(distance, taux_mut_bp, Tableau, taux_rec_bp):

    Couleurs = ['gold', 'coral', 'deepskyblue', 'limegreen']
    Styles = ['solid', 'dashed', 'dashdot', 'dotted' ]

    #Values that need to be changed if the table is not Tglobal
    nb_R = 4
    nb_pin = 2
    nb_f = 4

    #Retrieving the values of r to create the x-axis
    Echelle_R = Tableau[0:nb_R, 0]
    Echelle_distance = Echelle_R/taux_rec_bp

    #Computing the probabilities for each line
    les_probas, F, Pin = calcul_proba_nouvelle_mut(taux_mut_bp, Tableau, distance)

    #The case f=0 is treated appart from the others, as the value of p_in does not matter. We chose here pin = 0.5 (but this can be changed to take the set of datas given by pin=0.9
    plt.figure(figsize=(18,10))
    plt.plot(Echelle_distance, les_probas[0:nb_R], color = Couleurs[0], linestyle = Styles[0], linewidth=2, label = r'$f=%.d$' %(F[0]))
    #Then all the others curves are plotted
    for i in range(1,nb_f):
        plt.plot(Echelle_distance, les_probas[i*nb_R:(i+1)*nb_R], color = Couleurs[i], linestyle = Styles[0], linewidth=2, label = r'$f=%.1f, p_{in}=%.1f$' %(F[i*nb_R], Pin[0]))
        plt.plot(Echelle_distance, les_probas[nb_f*nb_R+i*nb_R:nb_f*nb_R+(i+1)*nb_R], color = Couleurs[i], linestyle = Styles[1], linewidth=2, label = r'$f=%.1f, p_{in}=%.1f$' %(F[i*nb_R+nb_f*nb_R], Pin[nb_f*nb_R]))


    #As r spans several order of magnitude, a log scale is chosen
    plt.xscale('log')

    plt.legend(prop={'size':15}, loc='lower left', ncol=1)

    plt.ylabel(r'Probability that a new mutation arises before purging', fontsize=20)
    plt.xlabel(r'Distance between the mating-type locus and the first mutation (bp)', fontsize=20)

    plt.show()
    plt.savefig("Fig5_ProbaNewMut.png")


#===========================
#Parameters and Commands
#===========================

#Note that Tglobal was imported before (line 27)

taux_mut_bp = 10**(-8)
taux_rec_bp = 10**(-6)
distance = 10**6
trace_probanewmut(distance, taux_mut_bp, Tglobal, taux_rec_bp)