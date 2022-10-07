import numpy as np

#This file to simulate the branching process and generate the datas needed to plot histograms

#============================
#PARTIAL DOMINANCE SCENARIO
#============================

#============================
#Functions to compute matrices used to implement the branching process, in the partial dominance case
#============================

def auxiliaire_a(r, Tet):
    return 1 - r*(0.75*Tet[0] + Tet[1]) + 0.25*Tet[1]*r**2

def auxiliaire_b(r, Tet):
    return r*(0.25*Tet[0] + 0.5*Tet[1] - 0.25*r*Tet[1])

def auxiliaire_c(r, Tet):
    return 0.25*r*(Tet[0] + r*Tet[1])


def matrice_A(R, Tet, S, f):

    A = np.zeros((3,3))

    A[0][0] =  S[0]*(f*auxiliaire_a(R[0], Tet) + (1-f)*(1-0.5*R[0]))
    A[0][1] =  S[0]*(f*auxiliaire_c(R[1], Tet) + (1-f)*(0.5*R[1]))
    A[0][2] =  S[0]*(1-f)

    A[1][0] =  S[1]*(f*auxiliaire_c(R[0], Tet) + (1-f)*(0.5*R[0]))
    A[1][1] =  S[1]*(f*auxiliaire_a(R[1], Tet) + (1-f)*(1-0.5*R[1]))
    A[1][2] =  S[1]*(1-f)

    A[2][0] =  S[2]*(f*auxiliaire_b(R[0], Tet))
    A[2][1] =  S[2]*(f*auxiliaire_b(R[1], Tet))
    A[2][2] =  S[2]*f

    return A

def matrice_T(R, Tet, S, f):

    T = np.zeros((3,3))

    T[0][0] = S[3]*(f*auxiliaire_b(R[0], Tet) + 0.5*(1-f))
    T[1][1] = S[3]*(f*auxiliaire_b(R[1], Tet) + 0.5*(1-f))

    return T

def matrice_D(R, Tet, S, f):

    return S[3]*np.eye(3)

#==============================
#Functions to compute the matrix of possible descendant vectors and the probability to produce each descendant vector
#==============================

#Matrice_enfants : Output : a matrix that displays, on each row i, the descendant vectors that an individual of type i can produce
#With three types of individuals (n=3), the matrix equals [[0, e1, e2, e3, e1+e1, e2+e1, e3+e1],[0, e1, e2, e3, e1+e2, e2+e2, e3+e2],[0, e1, e2, e3, e1+e3, e2+e3, e3+e3] ]

def Matrice_enfants(n):

    Enfants = [[0 for k in range(2*n+1)] for i in range(n)]

    #First column : null descendant vectors (no descendants)
    for k in range(n) :
        Enfants[k][0] = np.zeros(n)

    #Then a block with vectors e_i (one descendant)
    for k in range(n):
        for i in range(n):
            Enfants[k][i+1] = np.eye(n)[i]

    #Then a block with two descendants
    for k in range(n):
        for i in range(n):
            Enfants[k][n+1+i] = np.eye(n)[i] + np.eye(n)[k]

    return Enfants

#proba_enfants : Output : Matrix that contains the probability to produce each possible descendant vector (gives the probability to obtain each entry of Matrice_enfants)
#For each type i, descendants vectors are in the following order : [0, e1, e2, e3, e1+ei, e2+ei, e3+ei]

def proba_enfants(A, T, D):

    n=len(A[0])

    Tableau_probas = np.zeros((n, 2*n+1))

    # Filling the table with the rates
    for k in range(n):
        Tableau_probas[k,0] = D[k,k]

    for k in range(n):
        for i in range(n):
            Tableau_probas[k,1+i] = T[i][k]

    for k in range(n):
        for i in range(n):
            Tableau_probas[k,1+n+i] = A[i][k]

    #We normalize each entry by the sum of all rates to obtain a probability
    invTaux_individuel = 1/np.sum(A+T+D, axis = 0)
    invC = np.diag(invTaux_individuel) #Matrix that contains the inverted rates 1/c_i

    return invC.dot(Tableau_probas)

#===============================
#Branching process simulation
#===============================

#A function to have numerical simulations of the branching process, following a Gillespie algorithm
#R : list of the form [r1, r2, 1, 1] that specifies the recombination rates for each genotype. For the present paper, the same recombination rates were used for the two first genotypes
#Tet : list of the form [pin, 1-pin] that specifies the probability of intra-tetrad selfing
#s, h, f : selection coefficient, dominance coefficient, selfing rate
#X0 : list of the form [X01, X02, X03] that specifies the initial number of individuals of each genotype that carries the deleterious mutation.

#This function runs a single trajectory of the branching process, starting with X0 individuals of each genotype, and under the conditions specified by the parameters.
#The output is the time that it tool for the process to go extinct

def branchement(R, Tet, s, h, f, X0):

    #Initialization
    n = len(X0)
    X = X0
    Temps = 0

    #Computation of the matrices used to run the branching process, using auxiliary functions
    Survie = [1-h*s, 1-h*s, 1-s, 1]

    A = matrice_A(R, Tet, Survie, f)
    T = matrice_T(R, Tet, Survie, f)
    D = matrice_D(R, Tet, Survie, f)

    #Matrices used to update the composition of the population
    ME = Matrice_enfants(n)#Contains the possible descendant vectors
    PE = proba_enfants(A, T, D)#Contains the probability with which each descendant vector is chosen

    Taux_individuel = np.sum(A+T+D, axis = 0)

    #Loop for the evolution of the population
    while np.sum(X) != 0:

        #Total rate at which a reproduction event takes place
        Taux_total = np.sum(Taux_individuel*X)

        #Time to wait for the next reproduction event (exponential law of parameter Taux_total)
        temps_avant_saut_suivant = np.random.exponential(Taux_total)
        #Updating the time during which the process has survived
        Temps = Temps + temps_avant_saut_suivant

        #Choice of the type of parent that reproduces. type_parent is an index
        type_parent = np.random.choice(n, replace=True, p=Taux_individuel*X/np.sum(Taux_individuel*X))

        #Index of the descendant vector (vecteur_enfants) by which the parent is replaced
        c_parent= Taux_individuel[type_parent]
        vecteur_proba_enfants = PE[type_parent][:]
        indice_vecteur_enfants = np.random.choice(2*n+1,p = vecteur_proba_enfants)

        #Update of the population by retrieving one individual of the parent type, and adding the descendant vector that is found in the matrix ME
        X = X - np.eye(n)[type_parent,:] + ME[type_parent][indice_vecteur_enfants]

    return Temps


#==============================
#Function to lauch the branching process a certain number of times, with the same parameters
#==============================

def echantillon_temps_parametres(R, Tet, s, h, f, X0, nb_trajectoires):

    nb_parametres = 8
    les_Temps = np.zeros(nb_trajectoires+nb_parametres)

    #Saving the parameters at the beginning of the list that will contain the survival time of each run
    les_Temps[0] = R[0]#recombination rate
    les_Temps[1] = Tet[0]#pin
    les_Temps[2] = s
    les_Temps[3] = h
    les_Temps[4] = f
    les_Temps[5] = X0[0]#initial conditions (number of individuals of type 1)
    les_Temps[6] = X0[1]
    les_Temps[7] = nb_trajectoires


    for k in range(nb_trajectoires):

        les_Temps[k+nb_parametres] = branchement(R, Tet, s, h, f, X0)

    return les_Temps


#================================
#Commands used to generate datas for Figures 4 and S3
#================================

#Those commands will create tables containing 100 000 independant runs of the brancing process with the same parameters.
#ATTENTION : with 100 000 runs, the simulation can take several hours
#Figure 4
#r=0.001, pin=0.5, s=0.1, h=0.1, f=0.5, X0=[1,0,0], nb_trajectoires=100000
#T_r0001_s01 = echantillon_temps_parametres([0.001, 0.001, 1, 1], [0.5, 0.5], 0.1, 0.1, 0.5, [1, 0, 0], 100000)
#r=0.1, pin=0.5, s=0.1, h=0.1, f=0.5, X0=[1,0,0], nb_trajectoires=100000
#T_r01_s01 = echantillon_temps_parametres([0.1, 0.1, 1, 1], [0.5, 0.5], 0.1, 0.1, 0.5, [1, 0, 0], 100000)
#r=0.5, pin=0.5, s=0.1, h=0.1, f=0.5, X0=[1,0,0], nb_trajectoires=100000
#T_r05_s01 = echantillon_temps_parametres([0.5, 0.5, 1, 1], [0.5, 0.5], 0.1, 0.1, 0.5, [1, 0, 0], 100000)


#Figure S3
#r=0.001, pin=0.5, s=0.1, h=0.1, f=1, X0=[1,0,0], nb_trajectoires=100000
#T_f1_s01 = echantillon_temps_parametres([0.001, 0.001, 1, 1], [0.5, 0.5], 0.1, 0.1, 1, [1, 0, 0], 100000)
#r=0.1, pin=0.5, s=0.1, h=0.1, f=0.5, X0=[1,0,0], nb_trajectoires=100000
#T_f05_s01 = echantillon_temps_parametres([0.1, 0.1, 1, 1], [0.5, 0.5], 0.1, 0.1, 0.5, [1, 0, 0], 100000)
#r=0.5, pin=0.5, s=0.1, h=0.1, f=0, X0=[1,0,0], nb_trajectoires=100000
#T_f0_s01 = echantillon_temps_parametres([0.5, 0.5, 1, 1], [0.5, 0.5], 0.1, 0.1, 0, [1, 0, 0], 100000)

