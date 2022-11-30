# Instruction of LSCV-MCEA package
In this study, we proposed the concept of multi-objective optimization based structural network control principles (MONCP) by considering minimum driver nodes and maximum prior-known drug-target information.n. To solve MONCP, we formed a discrete multi-objective optimization problem with large-scale constrained variables and developed a
novel evolutionary optimization model (denoted LSCV-MCEA) by adapting multitasking framework and rankings-based fitness function method
This package includes Matlab scripts and several datasets for demo of LSCV_MCEA approach:

‘main_LSCV_MCEA.m’ is a Matlab function for the routine of experimental analysis. LSCV-MCEA aims to identify personalized drug targets by considering minimum driver
nodes and maximum prior-known drug-target information.

The input (case: BRCA) include:
(1)Path: The path of the user where the 'Main,m' is located.
(2)Cancer Data: tumor and normal sample data of BRCA.

The output results:
(1)PGIN_BRCA: BRCA patients' personalized gene interaction network construct by SSN method.
i.‘BRCA_i_PGIN.mat’ indicates that the PGIN of the i-th BRCA patient which contains the subnetwork adjacency matrix and the name of the gene in the subnetwork.
(2)BRCA_result: Non-dominated solutions of patient samples obtained by LSCV_MCEA.

i.‘BRCA_sample_i_LSCV_MCEA_boxchart.mat’ stores non-dominated solutions for the i-th patient by running the evolutionary algorithm 30 times each time.
ii.‘BRCA_sample_i_LSCV_MCEA_PF.mat’ indicates that a group Pareto front solutions of i-th patient obtained by performing non-dominated sort on boxchart.
iii.‘BRCA_sample_1_LSCV_MCEA_PS.mat’ indicates the drug targets corresponding to Pareto front solutions.
Suggestions

(1)Hardware suggestions for running this package: Window 10 or above; Matlab 2016 or above; RAM 32G or above.

(2)When users analyzed running this package, please note that:
i.Users should set the path in the program, firstly.
ii.Parameter setting of Popnum, Max_CalNum, and Experiment_num will affect the running time. With default parameters, LSCV_MCEA takes about 40 minutes to identify drug targets for a BRCA patient. Users can decrease running time by modifying above parameter.
iii.If users want to run their own experimental data, users should add the data of normal samples, tumor samples and mutation samples for the corresponding disease in the 'Code_construct_personalized_network' file.
iv.The user can modify the corresponding index in lines 44 and 57 of the 'DMOP_LSCV' function to identify one patient’s drug targets.

%    $Id: Main.m Created at 2022-11-28$ 
%   $Copyright (c) 2022 by School of Electrical Engineering, Zhengzhou University, Zhengzhou 450001, China$; 
%    $If any problem, please contact zero999396@outlook.com for help. $
