limit = [1, 10;1, 10]
for i = 1:2
    x(:,i) = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(100, 1);%初始种群的位置x
end

y=rand(10,5).*(10-5)+5;