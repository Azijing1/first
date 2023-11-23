clear all;
close all;
clc;
% 这份代码主要尝试了固定惯性值、线性递增惯性值、线性递减惯性值的区别
% f(x)=sum(x.*x-10*cos(2*pi.*x)+10)的最小值，极小点x=（0,0,...,0），理论上最小值分f(0,0,...,0)=0。
N=100;                         %群体粒子个数
D=3;                          %粒子维数
T=15;                         %最大迭代次数
c1=1.5;                        %学习因子1
c2=1.5;                        %学习因子2
w=0.8;                         %惯性权重
wmax=1;                      %惯性权重最大值
wmin=0.2;                      %惯性权重最小值
Xmax=100;                       %位置最大值
Xmin=-100;                      %位置最小值
Vmax=10000;                       %速度最大值
Vmin=-10000;                      %速度最小值
delta=0.00001;                    %事后估计法误差值
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
index=0;
%按照公式依次迭代直到满足精度或者迭代次数
for i=1:T
    %w=wmax-(wmax-wmin)*i/T; %%权值线性递减(最开始的算法)容易长时间陷入局部最优解
    %w=wmin+(wmax-wmin)*i/T; %%权值线性递增（最新的想法）
    %w=0;
    %psin=pi*(i-1)/(T-1)+pi/2;  w=((wmin-wmax)*sin(psin)+(wmin+wmax))/2;  %%非线性递增型
    wn(i)=w;
    for j=1:N
        if (func1(x(j,:))<pbest(j)) %%更新个体最优位置和个体最优值
            p(j,:)=x(j,:);
            pbest(j)=func1(x(j,:));
        end
        if (pbest(j)<gbest)         %%更新群体最优位置和群体最优值
            g=p(j,:);
            gbest=pbest(j);
        end
%         v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:));
%         x(j,:)=x(j,:)+v(j,:);
         x(j,:)=w*x(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:)); %%直接位置更新策略
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
    if gbest<15 && index==0
        index=i
    end
        
%     if i>1 && index==0
%         err=gb(i)-gb(i-1);
%         if err<delta
%             index=i
%         end
%     end   
end
g=g';                         %最优个体
gb(end);                   %最优值
figure(1)
plot(gb)
%ylim([0,80]);
xlabel('迭代次数')
ylabel('适应度值')
title('适应度进化曲线')
%适应度函数
% figure(2);
% ix=1:T;
% plot(ix,wn,ix,wn1);
% xlim([-2,230]);
% ylim([0.4,1.4]);
% legend('线性递增w','非线性递增w');
function result=func1(x)
%summ=sum(x.^2); %先对矩阵每个数进行平方，再对矩阵的每个数求和
summ=sum(x.*x-10*cos(2*pi.*x)+10);
%summ=sum(x);
result=summ;
end


