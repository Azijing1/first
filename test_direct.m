clear all;
close all;
clc;
%%直接位置型更新策略测试
n=40;
index=0;

for i=2:2
    for j=1:2
        a=zeros(1,n);
        b=zeros(1,n);
        for k=1:n
            a(k)=direct_partcle1(i,j);%%0线性递减型（左边）                              1选择速度型位置更新策略（右边）
            %b(k)=direct_partcle1(2,2);%%1线性递增型（左边）                              2或者其他选择直接位置型位置更新策略（右边）
        end                           %%2非线性递增型 ，其他数字w恒定为0.8（左边）
            index=index+1;
             rec_best(index,:)=a;
        %   xlim([-2,2]);
        %   ylim([-2,1200]);
            %mean(b)
        
    end
end
% figure(1);
% ix=1:n;
% plot(ix,rec_best(1,:),ix,rec_best(2,:),ix,rec_best(3,:),ix,rec_best(4,:));
% legend('线性递减型惯性系数最优值','线性递增型惯性系数最优值','非线性递增型惯性系数最优值','w恒定为0.8惯性系数最优值');
% title('速度位置型位置更新策略');
figure(1);
ix=1:n;
plot(ix,rec_best(1,:),ix,rec_best(2,:));
legend('速度-位置型下的最优值','位置直接型下的最优值');
title('非线性递增惯性权值下位置更新策略');

