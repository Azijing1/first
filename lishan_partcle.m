clc
clear
% 有N件物品和一个容量为W的背包。第i件物品的体积是volume[i]，价值是value[i]。求解将哪些物品装入背包可使价值总和最大同时又不能超过背包的总容积。
% 
% 其中：
% N = 10；
% W= 300；
% volume=[95 75 23 73 50 22 6 57 89 98]；
% value=[89 59 19 43 100 72 44 16 7 64]；

% 初始化种群
% IndN = 50; 
narvs = 10; %维数是十
n = 100; 
x_ub = 3;
x_lb = -3;

vmax = 1.2; %粒子的最大速度
w = 0.9;  % 惯性权重
c1 = 2;  % 每个粒子的个体学习因子，也称为个体加速常数
c2 = 2;  % 每个粒子的社会学习因子，也称为社会加速常数

%初始种群
x = randsrc(n,narvs,[0,1;0.5,0.5]);%产生一个n*narvs 100*10的矩阵，矩阵里的数字只有0和1，
                                   %且两者产生的概率各为0.5

v = -vmax + 2*vmax .* rand(n,narvs); %产生-1.2-1.2的随即速度二维表

%计算种群适应度
fitness = targetPackage(x',n);   %计算各个粒子的重量也就是适应度值
pbest = x;   % 初始化这n个粒子迄今为止找到的最佳位置（是一个n*narvs的向量）
ind = find(fitness == max(fitness), 1);  % 找到适应度最大的那个粒子的下标
gbest = x(ind,:);  % 定义所有粒子迄今为止找到的最佳位置（是一个1*narvs的向量）

K=100;  %迭代次数
fitnessbest=zeros(K,1); %保存每一代的函数值
best=0;
bestOne = zeros(1,10);

for t = 1:K    
    for i = 1:n
        v(i,:) = w*v(i,:) + c1*rand(1)*(pbest(i,:) - x(i,:)) + c2*rand(1)*(gbest - x(i,:));  %因为括号里是最优值减去当前值，所以可以近似认为速度很大时，最优点为1当前点为0
                                                                                             %速度很小可以近似认为，最优点为0，当前点为1，速度中等时两者同1或者同0
                                                                                             %这就能解释下面为何，速度大时要，让当前点靠近1，速度小时，让当前点靠近0
        % 判断速度是否超过限制了(还可以取模)
        for j = 1:narvs
            if v(i,j) < -vmax
                v(i,j) = -vmax;
            elseif v(i,j) > vmax
                v(i,j) = vmax;
            end
        end
        %sigmoid函数将例子的速度映射到0-1之间
        vs(i,:)=1./(1+exp(-v(i,:))); 
        for j = 1:narvs    %二元离散的关键步骤
            if rand < vs(i,j)  %奇怪的问题，上1下0从多次结果来看收敛更快，也基本可以收敛到最好388左右
                x(i,j) = 1;    %但上0下1也有较大机会收敛出388，偶尔368左右，但收敛很慢（已解决）
            else
                x(i,j) = 0;
            end
        end
        
        fit = targetPackage(x(i,:)',1);  % 重新计算第i个粒子的适应度
        if fit > targetPackage(pbest(i,:)',1)   % 如果第i个粒子的适应度大于这个粒子迄今为止找到的最佳位置对应的适应度
            pbest(i,:) = x(i,:);   % 那就更新第i个粒子迄今为止找到的最佳位置
        end
        if  fit > targetPackage(gbest',1)  % 如果第i个粒子的适应度大于所有的粒子迄今为止找到的最佳位置对应的适应度
            gbest = pbest(i,:);   % 那就更新所有粒子迄今为止找到的最佳位置
        end
    end
    h = targetPackage(gbest',1);    
    if h>best         
        best = h;  
        bestOne = gbest; 
    end
    fitnessbest(t,1)=best;
end
bestOne
plot(1:K,fitnessbest,'-');
grid on;
grid minor;
