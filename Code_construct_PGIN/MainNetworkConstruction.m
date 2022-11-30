clc
clear
%************************part1:LOAD sample data and network data************************
%*********Case 1:BRCA data**********
unzip('BRCA_normal.zip') 
unzip('BRCA_tumor.zip') 
expression_tumor_fileName = 'BRCA_tumor.txt';
expression_normal_fileName = 'BRCA_normal.txt';
unzip('BRCA_mutation_data.zip') 
load('BRCA_mutation_data.mat')

% %*********Case 2:LUSC(1-49)+LUAD data(50-106)**********
% unzip('LUNG_normal.zip') 
% unzip('LUNG_tumor.zip') 
% expression_tumor_fileName = 'LUNG_tumor.txt';
% expression_normal_fileName = 'LUNG_normal.txt';
% unzip('LUNG_mutation_data.zip') 
% load('LUNG_mutation_data.mat')


% *************************tumor data****************************

[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);
Sample_name_tumor=tumor.textdata(1,2:end);
Tumor=tumor.data;

%*************************normal data****************************

[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);
Normal=normal.data;

%********************part2:SSN algorithm ***************************************
%read the size of data
new_T=Tumor;
new_N=Normal;
new_gene=gene_list;

[row,colunm]=size(new_T);
Ref=new_N; 

for i=1:size(new_T,2)
% for i=1:size(new_T,2)
    i
    tic     

    %for the i-th sample
    %construct the tumor SSN
    
    sample_tumor=new_T(:,i);
    [R0,P]=SSN(sample_tumor,Ref);
    
    P(isnan(P))=0;
    P(abs(P)>=0.05)=0;
    P(P~=0)=1;

    %construct the normal SSN
    clear  sample_tumor 
    sample_normal=new_N(:,i);
    [R1,P1]=SSN(sample_normal,Ref);
    clear  sample_normal 
    
    P1(isnan(P1))=0;
    P1(abs(P1)>=0.05)=0;
    P1(P1~=0)=1;
    
    C=P-P1; 
    R=abs(log2(abs(R0./R1)));
    R(isnan(R))=0;
    
    %sample_network{i,1}=C;
    D=abs(C).*R;
   %*****************save network data**************************8

   CC=D.*New_A;
   CC(isnan(CC))=0;CC(CC==inf)=0;

   k=sum(CC);%
   [b,a]=find(k~=0);
    
   subnetwork_nodes=a';
   subnetwork_genes=gene_list(subnetwork_nodes);

   scores=k(a)';

   subnetwork_adjacency0=CC(subnetwork_nodes,:);
   subnetwork_adjacency=subnetwork_adjacency0(:,subnetwork_nodes);

   subnetwork_adjacency(subnetwork_adjacency~=0)=1;
   name=['Personalized_sample_data_' expression_tumor_fileName(1:4) '_'];
   samplename=strcat(name,num2str(i),'.mat');
   save(samplename,'subnetwork_genes','subnetwork_adjacency')

    
    
    toc 
    
end