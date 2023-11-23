clear all;
close all;
clc;
% 例1：计算函数
% f(x)=∑（i=1：n）xi^2(−20≤xi≤20)的最小值，其中个体x的维数n=10，这是一个简单的平方和函数，只有一个极小点x=（0,0,...,0），理论上最小值分f(0,0,...,0)=0。
N=100;                         %群体粒子个数
D=10;                          %粒子维数
T=200;                         %最大迭代次数
c1=1.5;                        %学习因子1
c2=1.5;                        %学习因子2
w=0.8;                         %惯性权重
Xmax=20;                       %位置最大值
Xmin=-20;                      %位置最小值
Vmax=10;                       %速度最大值
Vmin=-10;                      %速度最小值
%初始化个体
x=rand(N,D)*(Xmax-Xmin)+Xmin; %横坐标一排是x1-x10,纵坐标一列就是1-100粒子，每个粒子的一行值可计算出一个可能是最优的值
v=rand(N,D)*(Vmax-Vmin)+Vmin; %产生一系列Vmin-Vmax的随机数
%初始化个体最优位置和最优值
p=x;
pbest=ones(N,1);
for i=1:N
    pbest(i)=func1(x(i,:)); %先对矩阵每个数进行平方，再对矩阵的每个数求和，
                           %其实这个就是每一行都算出一个平方求和值，一共一百个这样的值，就是对于一百种不同的x1-x10，求一百次f值
end
%初始化全局最优位置和最优值
g=ones(1,D);
gbest=inf;
for i=1:N
    if (pbest(i)<gbest)
        g=p(i,:);
        gbest=pbest(i);
    end
end
gb=ones(1,T);
%按照公式依次迭代直到满足精度或者迭代次数
for i=1:T
    for j=1:N
        if (func1(x(j,:))<pbest(j)) %%更新个体最优位置和个体最优值
            p(j,:)=x(j,:);
            pbest(j)=func1(x(j,:));
        end
        if (pbest(j)<gbest)         %%更新群体最优位置和群体最优值
            g=p(j,:);
            gbest=pbest(j);
        end
        v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:));
        x(j,:)=x(j,:)+v(j,:);
        %边界条件处理
        for ii=1:D
            if (v(j,ii)<Vmin)||(v(j,ii)>Vmax)
                v(j,ii)=rand*(Vmax-Vmin)+Vmin;
            end
            if (x(j,ii)<Xmin)|(x(j,ii)>Xmax)
                x(j,ii)=rand*(Xmax-Xmin)+Xmin;
            end
        end
    end
    %记录全局最优值
    gb(i)=gbest;
end
g;                         %最优个体
gb(end);                   %最优值
figure
plot(gb)
xlabel('迭代次数')
ylabel('适应度值')
title('适应度进化曲线')
%适应度函数
function result=func1(x)
summ=sum(x.^2); %先对矩阵每个数进行平方，再对矩阵的每个数求和
result=summ;
end


