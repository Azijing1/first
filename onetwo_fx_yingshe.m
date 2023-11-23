clear all;
close all;
clc;
% 例3：用离散粒子群算法求函数 f(x)=x+6sin(4x)+9cos(5x) 的最小值，其中x的取值范围是[0,9]，这是一个有局部多个极值的函数，图形如图所示:
N=100;                         %群体粒子个数
D=20;                           %粒子维数
T=200;                         %最大迭代次数
c1=1.5;                        %学习因子1
c2=1.5;                        %学习因子2
Wmax=0.8;                      %惯性权重最大值
Wmin=0.4;                      %惯性权重最小值
Xs=9;                        %位置最大值
Xx=0;                       %位置最小值
Vmax=10;                        %速度最大值
Vmin=-10;                       %速度最小值
global x_r
%初始化个体
x=randi([0,1],N,D);
v=rand(N,D)*(Vmax-Vmin)+Vmin;
%初始化个体最优位置和最优值
p=x;
pbest=ones(N,1);
for i=1:N
    pbest(i)=func3(x(i,:),Xs,Xx);
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
        %更新个体最优位置和最优值
        if (func3(x(j,:),Xs,Xx)<pbest(j))
            p(j,:)=x(j,:);
            pbest(j)=func3(x(j,:),Xs,Xx);
        end
           %更新全局最优位置和最优值
        if (pbest(j)<gbest)
            g=p(j,:);
            gbest=pbest(j);
            x_best=x_r
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
m=0;
gb(end)
for j=1:D
    m=g(j)*2^(j-1)+m;
end
f1=Xx+m*(Xs-Xx)/(2^D-1);%最优值
figure
plot(gb)
xlabel('迭代次数')
ylabel('适应度值')
title('适应度进化曲线')
%适应度函数
function result=func3(x,Xs,Xx)
global x_r
m=0;
D=length(x);
for j=1:D
    m=x(j)*2^(j-1)+m;
end
f=Xx+m*(Xs-Xx)/(2^D-1);%译码成十进制数,到这里可以将一行01矩阵转化为一个0-9的随机数
x_r=f;
fit=f+6*sin(4*f)+9*cos(5*f);
result=fit;
end
