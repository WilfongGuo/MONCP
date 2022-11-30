function Fitness=Ranking_Fitness(CV,f2,gen,VAR,N1,maxGen,HZ)


if HZ==1
    F2=f2(:,1);
else
    F2=f2(:,2);
end

feasiIndex=find(CV==0);
for i = 1:size(f2,1)
    if CV(i) <= VAR
        D11 = find((CV<=VAR)&(F2<=F2(i)));
        Rc1(i) = size(D11,1);
    else D21 = find((CV<=CV(i)));
        Rc1(i) = size(D21,1);
    end
end
%% 按照目标函数
for  k = 1:size(f2,1)
    E1 = find(F2<=F2(k));
    Rf1(k) = size(E1,1);
end
%% 选择
Pfea1(gen) = size(feasiIndex,1)/N1;

a(gen)=0.5+0.5*(Pfea1(gen)+gen/maxGen);
if a(gen)>1
    a(gen)=1;
end
Fitness = a(gen)*Rc1+(1-a(gen))*Rf1;
end


