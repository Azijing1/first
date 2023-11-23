clear all;
close all;
clc;
% 例4：0-1背包问题，有N件物品和容积为V的包，第i件物品的容积是c(i),价值是w(i)，求将这些物品放入包中，使物体的总容积不超过背包容积，且总价值和最大。
% 假设物品数量为10，背包容量为300.每件物品的体积为[95,75,23,73,50,22,6,57,89,98];价值为[89,59,19,43,100,72,44,16,7,64];
N=100;                          %群体粒子个数
D=10;                           %粒子维数
T=200;                          %最大迭代次数
c1=1.5;                         %学习因子1
c2=1.5;                         %学习因子2
Wmax=0.8;                       %惯性权重最大值
Wmin=0.4;                       %惯性权重最小值
Vmax=10;                        %速度最大值
Vmin=-10;                       %速度最小值
V=300;                          %背包容量
C=[95,75,23,73,50,22,6,57,89,98];            %物品体积
W=[89,59,19,43,100,72,44,16,7,64];           %物品价值
afa=2;                                       %惩罚函数系数
%初始化个体
x=randi([0,1],N,D);
v=rand(N,D)*(Vmax-Vmin)+Vmin;
%初始化个体最优位置和最优值
p=x;
pbest=ones(N,1);
for i=1:N
    pbest(i)=func4(x(i,:),C,W,V,afa);
end
%初始化全局最优位置和最优值
g=ones(1,D);
gbest=eps;
for i=1:N
    if (pbest(i)>gbest)
        g=p(i,:);
        gbest=pbest(i);
    end
end
gb=ones(1,T);
%按照公式依次迭代直到满足精度或者迭代次数
for i=1:T
    for j=1:N
         %更新个体最优位置和最优值
        if (func4(x(j,:),C,W,V,afa)>pbest(j))
            p(j,:)=x(j,:);
            pbest(j)=func4(x(j,:),C,W,V,afa);
        end
        %更新全局最优位置和最优值
        if (pbest(j)>gbest)
            g=p(j,:);
            gbest=pbest(j);
        end
         %计算动态惯性权重值
          w=Wmax-(Wmax-Wmin)*i/T;
          %更新位置和速度
          v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:));
     %边界条件处理
        for ii=1:D
             if (v(j,ii)<Vmin)||(v(j,ii)>Vmax)
             v(j,ii)=rand*(Vmax-Vmin)+Vmin;
             end
        end
        vx(j,:)=1./(1+exp(-v(j,:)));
        for jj=1:D
             if vx(j,jj)>rand
             x(j,jj)=1;
             else
             x(j,jj)=0;
             end
        end
    end
    gb(i)=gbest;
end
g;%最优个体
figure
plot(gb)
xlabel('迭代次数')
ylabel('适应度值')
title('适应度变化曲线')
%适应度函数
function result=func4(f,C,W,V,afa)
fit=sum(f.*W);
TotalSize=sum(f.*C);
if TotalSize<=V
    fit=fit;
else
    fit=fit-afa*(TotalSize-V);
end
result=fit;
end