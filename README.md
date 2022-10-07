# Mutations_near_mating_type
Code and datas for figures of the article "The fate of deleterious or overdominant mutations near mating-type loci under partial selfing"


The .png files contain the figures displayed in the paper.
The .py files contain the codes used to produce data (simulate the branching process) and plot figures
The files named 'Thisto....' or 'Tableau_concatene" are text files that contain the realisation of a certain number of independant run of the branching process, for a chosen selection scenario and a chosen parameter set (see the parameters used in the .py files named 'SimuBranchingProcess...'). 
The tables 'Thisto...' contain one row, with the first eigth entries containing the parameters used, and the rest of the row the values of extinction time for each run. 
'Tableau_concatene' contains multiple rows, with the first eight entries being the parameters used, and the rest of the rows the number of reproduction events that occurred for each run.

Figure 2 :
Python file containing codes and commands to produce the figure : Figure2_Heatmaps_VapOverdom.py
Saved figure : Fig2_Heatmaps_Overdom.png

Figure 3 : 
Python file containing codes and commands to produce the figure : Figure3_DeriveeVap.py
Saved figure : Fig3_derivee_contresel_overdom.png

Figure 4 (and S3):
Python file containing codes and commands to produce the figure : Figure4_S3_BarChart_PartialDominance.py
Data files needed to execute the program : Thisto_r0001_s01, Thisto_r01_s01, Thisto_r05_s01, Thisto_f0_s01, Thisto_f1_s01, Thisto_f05_s01
Generated with file : SimuBranchingProcess_PartialDominance.py
Saved figure : Fig4_Histo_contresel_r_s01.png and FigureS3_Histo_contresel_f.png

Figure 5 :
Python file containing codes and commands to produce the figure : Figure5_ProbaNewMut.py
Data files nedded to execute the program : Tableau_concatene
Produced with file : SimuBranchingProcess_ProbaNewMut_Fig5
Saved figure : Fig5_ProbaNewMut.png

Figure S4 (and S5):
Python file containing codes and commands to produce the figure : FigureS4_S5_BarChart_Overdominance.py
Data files needed to execute the program : T_histo_overdom_r01_s40001, T_histo_overdom_r05_s4001, T_histo_overdom_r05_s40001, Thisto_overdom_f01, Thisto_overdom_f05, Thisto_overdom_f09, Thisto_overdom_f1
Generated with file : SimuBranchingProcess_Overdominance.py
Saved figure : FigS4_Histo_overdom_r.png and FigS5_Histo_overdom_f.png
