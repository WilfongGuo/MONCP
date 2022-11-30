function [PF,PS,Non_dominated_sol]=LSCV_MCEA_EA(adjacent,N,MaxFE,experiment_num,label)
%%     input:
%                   adjacent            :    patient sample's PGIN
%                   N          :    population size
%                   MaxFE      :    the maximum number of function evaluation
%                   experiment_num  :    the number of algorithm rus
%       output:
%                   PF   :    the value of objective function of DNB obtained by MMPDNB
%                   PS   :    non dominated DNB obtained by MMPDNB
%                   NDS  :    the results of 30 runs of MMPDNB
%% MDS
Cons=adjacent.subnetwork_adjacency;
D=length(label); %维度
Cons(all(Cons==0,2),:)=[];%删除全零行
Cons=unique(Cons,'rows');

%% DFVS
% A=adjacent.subnetwork_adjacency;
% [z2,z1]=find(A~=0);
% z=[z1 z2];
% N1=length(A);
% [N2,~]=size(z);
% %calculate the adjacency matrix of bipartite graph
% A_adjacent=zeros(N2,2*N1);
% for i=1:N2
%     
%     A_adjacent(i,z(i,1))=1;
%     A_adjacent(i,z(i,2))=-1;
%     A_adjacent(i,N1+z(i,1))=N1;
%     
% end
% Cons=A_adjacent;
% Cons(all(Cons==0,2),:)=[];%删除全零行
% Cons=unique(Cons,'rows');

%% NCUA
% Network=adjacent.subnetwork_adjacency;
% [z1,z2]=find(Network~=0);
% z=[z1 z2];
% N1=length(Network);
% [N2,~]=size(z);%calculate the adjacency matrix of bipartite graph
% A_adjacent=zeros(N2,N1);
% for i=1:N2
%     
%     A_adjacent(i,z(i,1))=1;
%     A_adjacent(i,z(i,2))=1;
% end
% Cons=A_adjacent;
% Cons(all(Cons==0,2),:)=[];%删除全零行
% Cons=unique(Cons,'rows');


N1=N*0.3;    % the population size of subpopulations
maxGen=2*ceil(MaxFE/(N+2*N1));
%% Set experiment parameters
Non_dominated_sol=cell(experiment_num,1);

for EXP_NUM=1:experiment_num
    
    X=0;
    gen=1;
    %% Generate intial population.
    FE=0;   % the number of function evaluation
    Population1=Initialization(D,N);
    Population2=Initialization(D,N1);
    Population3=Population2;
    
    f1=Calfunctionvalue(Population1,label);% Calculate the objective function of 
    CV1=Calcons(Population1,Cons);% Calculate the Constraint violation degree   
    f2=Calfunctionvalue(Population2,label);
    CV2=Calcons(Population2,Cons);   
    f3=Calfunctionvalue(Population3,label);
    CV3=Calcons(Population3,Cons);
    
    VAR0=max([CV1;CV2;CV3]);    
    Fitness1 = CalFitness(f1,CV1,VAR0);%Calculate the Fitness
    Fitness2=Ranking_Fitness(CV2,f2,gen,VAR0,N1,maxGen,1);
    Fitness3=Ranking_Fitness(CV3,f3,gen,VAR0,N1,maxGen,2);   
    FE=FE+N+2*N1;
      
    %% Evolution
    while(FE<MaxFE)
               
        cp=(-log(VAR0)-6)/log(1-0.5);
        %             adjust the threshold
        if X < 0.5
            VAR=VAR0*(1-X)^cp;
        else
            VAR=0;
        end
        %% Tournament Selection
        MatingPool1= TournamentSelection(2,N,Fitness1);
        MatingPool2= TournamentSelection(2,N1,Fitness2);
        MatingPool3= TournamentSelection(2,N1,Fitness3);
        %% Generate Offspring
        Offspring1=OperatorGAhalf(Population1(MatingPool1));
        Offspring2=OperatorGAhalf(Population2(MatingPool2));
        Offspring3=OperatorGAhalf(Population3(MatingPool3));                            
        %% Tournament Selection
        pop_new1=[Population1;Offspring1;Offspring2;Offspring3];
        f_new=Calfunctionvalue(pop_new1,label);
        CV_new=Calcons(pop_new1,Cons);
        Fitness_new1 = CalFitness(f_new,CV_new,VAR);
        FE=FE+size(Offspring1,1);
        [Population1,Fitness1] = EnvironmentalSelection(pop_new1,N,Fitness_new1); % population1
                     
        pop_new2=[Population2;Offspring1;Offspring2;Offspring3];
        f_new2=Calfunctionvalue(pop_new2,label);
        CV_new2=Calcons(pop_new2,Cons);
        Fitness_new2=Ranking_Fitness(CV_new2,f_new2,gen,VAR0,N1,maxGen,1);
        FE=FE+size(Offspring2,1);
        [~,rank]   = sort(Fitness_new2);
        Population2 = pop_new2(rank(1:N1));
        Fitness2=Fitness_new2(rank(1:N1));% population2
   
        pop_new3=[Population3;Offspring1;Offspring2;Offspring3];
        f_new3=Calfunctionvalue(pop_new3,label);
        CV_new3=Calcons(pop_new3,Cons);
        Fitness_new3=Ranking_Fitness(CV_new3,f_new3,gen,VAR0,N1,maxGen,2);
        FE=FE+size(Offspring3,1);
        [~,rank]   = sort(Fitness_new3);
        Population3 = pop_new3(rank(1:N1));
        Fitness3=Fitness_new3(rank(1:N1));% population3
        
        X=X+1/maxGen;
        gen=gen+1;
    end   
    %% Record non-dominated soutions
    f1=Calfunctionvalue(Population1,label);
    CV1=Calcons(Population1,Cons);
    ind=find(CV1==0);
    Population1=Population1(ind,:);
    objs=f1(ind,:);   
    [FrontNo,~] = NDSort(objs,inf);
    outputpop=Population1(FrontNo==1,:,:);
    Non_dominated_sol{EXP_NUM}=outputpop;
end

Non_dominated_sol=cell2mat(Non_dominated_sol);

[pop,~,~]=unique(Non_dominated_sol,'rows');% De-duplication

functionvalue = Calfunctionvalue(pop,label);

[FrontNo,~] = M_non_domination_scd_sort(pop,functionvalue);

POP=pop(FrontNo==1,:);

FV=functionvalue(FrontNo==1,:);

[PF,mod_position]=sortrows(FV);

PS=POP(mod_position,:);
end

function CV=Calcons(pop,Cons)
pp=size(pop,1);
    %%  约束
    CV=zeros(pp,1);
    for z=1:pp
        ind=pop(z,:).';
        cv=Cons*ind;
        index=find(cv==0);
        CV(z)=length(index);
    end
end

function pop =Initialization(D,N)
pop=rand(N,D);
for i = 1 : N   
    for j = 1 : D
        if pop(i,j)<=0.5
            pop(i,j)=0;
        else
            pop(i,j)=1;
        end
    end   
end
end
function f=Calfunctionvalue(pop,label)
%% 函数值
pp=size(pop,1);
f=zeros(pp,2);
    for i= 1:pp
        precision(i,:)=pop(i,:)*label;
        x(i,:)=sum(pop(i,:),2);
    end
    
    f(:,1)=x/length(label);
    f(:,2)=1-precision/length(label);
end

function index = TournamentSelection(K,N,varargin)
    varargin    = cellfun(@(S)reshape(S,[],1),varargin,'UniformOutput',false);
    [Fit,~,Loc] = unique([varargin{:}],'rows');
    [~,rank]    = sortrows(Fit);
    [~,rank]    = sort(rank);
    Parents     = randi(length(varargin{1}),K,N);
    [~,best]    = min(rank(Loc(Parents)),[],1);
    index       = Parents(best+(0:N-1)*K);
end

function Offspring = OperatorGAhalf(Parent)
 proM=1;
Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
            %% Genetic operators for binary encoding
            % Uniform crossover
            k = rand(N,D) < 0.5;

            Offspring1    = Parent1;
            Offspring2    = Parent2;
            Offspring1(k) = Parent2(k);
            Offspring2(k) = Parent1(k);
            Offspring     = [Offspring1;Offspring2];
            if size(Offspring,1)<size(Parent,1)
                ofd=Parent1(3,:);
                ofd(1,k(3,:))=Parent2(3,k(3,:));
                Offspring= [Offspring;ofd];
            end
            % Bit-flip mutation
            Site = rand(size(Offspring,1),D) < proM/D;
            Offspring(Site) = ~Offspring(Site);
end

function [Population,Fitness] = EnvironmentalSelection(Population,N,Fitness)
    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        pop=Population(Next,:);
        objs=Calfunctionvalue(pop,label);
        Del  = Truncation(objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next,:);
    Fitness    = Fitness(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank,:);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end

