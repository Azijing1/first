function fitness = targetPackage(x,indNum)
    
    volume=[95 75 23 73 50 22 6 57 89 98];        %物品体积
    value=[89 59 19 43 100 72 44 16 7 64];          %物品价值
    Weight=300;  %背包重量
    a = zeros(indNum,1);
    for i=1:indNum
        a(i,1) = volume*x(:,i); %计算每个粒子（）100个粒子的体积
    end
    %超过总重量的个体的适应度值都视为0
    k = find(volume*x<Weight); %查找哪些粒子的重量没有超标并给出他们的索引值
    fitness=zeros(indNum,1);
    for j=1:size(k,2)  %size(k,2)表示k的列数
        fitness(k(j),1) = value*x(:,k(j));   %记录没有超重粒子的重量，其余为0
    end   
end
