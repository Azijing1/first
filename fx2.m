clear all;
close all;
clc;
x=-4:0.01:4;
y=-4:0.01:4;
N=size(x,2);     %输出x第二维度的大小 这里就是列数
for i=1:N
    for j=1:N
        z(i,j)=3*cos(x(i)*y(j))+x(i)+y(j)*y(j);
    end
end
mesh(x,y,z)
xlabel('x')
ylabel('y')
