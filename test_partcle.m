clear all;
close all;
clc;
n=40;
a=zeros(1,n);
b=zeros(1,n);
for i=1:n
    a(i)=partcle1(3);%%0线性递减型
    b(i)=partcle1(2);%%1线性递增型
end                  %%2非线性递增型 ，其他数字w恒定为0.8
    figure(1);
    ix=1:n;
    plot(ix,a,ix,b);
%     xlim([-2,23]);
%     ylim([0,150]);
    legend('恒定惯性系数最优值','非线性递增惯性系数最优值');
    mean(a)
    mean(b)