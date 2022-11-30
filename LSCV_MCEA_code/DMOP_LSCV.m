function gen_name=DMOP_LSCV(expression_tumor_fileName,expression_normal_fileName,path,N,MaxFE,experiment_num)
%% *************************tumor data****************************
[tumor,~,~]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);
Tumor=tumor.data;

%% *************************normal data***************************
[normal,~,~]=importdata(expression_normal_fileName);
Normal=normal.data;

%% *************************network data**************************
if isequal(expression_tumor_fileName(1:4),'BRCA')
    Z=load([path,'Code_construct_personalized_network\',expression_tumor_fileName(1:4),'_mutation_data.mat']);
    Z=Z.Z;
elseif isequal(expression_tumor_fileName(1:4),'LUNG')
    Z=load([path,'Code_construct_personalized_network\',expression_tumor_fileName(1:4),'_mutation_data.mat']);
    Z=Z.Z;
end

N1=length(gene_list);
[N2,~]=size(Z);
Net=zeros(N1);
for i=1:N2
    
        Net(Z(i,2),Z(i,1))=1;
        Net(Z(i,1),Z(i,2))=1;
end

%calculate the adjacency matrix of PPI 

N1=length(gene_list);
[N2,~]=size(Z);
Net=zeros(N1);
for i=1:N2
    
        Net(Z(i,2),Z(i,1))=1;
        Net(Z(i,1),Z(i,2))=1;
end
%% Construct PGIN and save 
path1=strcat(path,'PGIN_',expression_tumor_fileName(1:4),'\');
mkdir(path1);
Ref=Normal; 

for i=1:size(Tumor,2) % all patient samples
    [subnetwork_genes,subnetwork_adjacency] = construct_PGIN(i,Normal, Tumor,gene_list,Net,Ref);
    
    filename=strcat(path1,expression_tumor_fileName(1:4),'_',num2str(i),'_PGIN.mat');
    save(filename,'subnetwork_genes','subnetwork_adjacency');
end


%% EA
path2=strcat(path,expression_tumor_fileName(1:4),'_result','\');
mkdir(path2); % output path

numf = dir(path1);
for fnum = 1:length(numf)-2 % all patient samples
    adjacent=load([path1,expression_tumor_fileName(1:4),'_',num2str(fnum) '_PGIN']);
    unzip('label.zip');              
    b=load(strcat('BRCA_',num2str(fnum)));% label
    label=b.m;
    [PF,PS,Non_dominated_sol]=LSCV_MCEA_EA(adjacent,N,MaxFE,experiment_num,label);
    
    filename=strcat(path2,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_LSCV_MCEA_PF.mat');
    save(filename,'PF');
    filename2=strcat(path2,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_LSCV_MCEA_PS.mat');
    save(filename2,'PS');
    filename3=strcat(path2,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_LSCV_MCEA_boxchart.mat');
    save(filename3,'Non_dominated_sol');
end
%% output 

for fnum=1:length(numf)-2 % all patient samples
    data=load([path1,expression_tumor_fileName(1:4),'_',num2str(fnum) '_PGIN']);
%     Sol=load([path2,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_LSCV_MCEA_PF.mat']);
%     Sol=Sol.PF;
    POP=load([path2,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_LSCV_MCEA_PS.mat']);
    POP=POP.PS;
    

    
    LSCV_MCEA=cell(size(POP,1),1);
    for i=1:size(POP,1)
        P=find(POP(i,:)~=0);
        for j=1:length(P)
            MM{j,1}=data.subnetwork_genes{P(j),1};
        end
        LSCV_MCEA{i,1}=MM;
        clear MM
    end
    gen_name{fnum,1}=LSCV_MCEA;
end
end