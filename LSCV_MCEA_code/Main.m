clc
clear

                 
%% ***************************Input*************************************

path='...\';        % The path of 'Main,m' on the user's computer and '\' need reserve.
%BRCA data
unzip('BRCA_tumor.zip');              % Take BRCA as an example
unzip('BRCA_normal.zip');
 
expression_tumor_fileName = strcat('BRCA_tumor.txt');
expression_normal_fileName = strcat('BRCA_normal.txt');

%Another Lung data
%unzip('LUNG_tumor.zip');              % Take LUNG as an example
%unzip('LUNG_normal.zip');
 
%expression_tumor_fileName = strcat('LUNG_tumor.txt');
%expression_normal_fileName = strcat('LUNG_normal.txt');


%% *******************The function of DMOP_LSCV****************************

%% default parameters
popnum=300;
Max_CalNum=100000;
Experiment_num=30;

gen_name=DMOP_LSCV(expression_tumor_fileName,expression_normal_fileName,path,popnum,Max_CalNum,Experiment_num);

