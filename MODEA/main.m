%% 设置多次实验
clc 
clear
close all
num_experiments = 1;
%% 算法数据
NP = 200;          % 种群数量
maxgen_base = 200;     % 迭代次数
Pc = 0.8;
Pm = 0.2;
M = 2;            % 目标函数个数
rho = 0.5; % 最大收益目标函数买家花费占比
ka = 10;  
v = 0.5;    % 最大收益目标函数卖家收入占比
epsilon = 1;
sigma = 174;
omiga = 40; %带宽
p = 100; % 传输功率
T_delay = zeros(num_experiments,1); % 总延迟时间
delay_average_results = zeros(num_experiments,1); % 平均延迟时间
x = zeros(num_experiments,1); % 车辆数横坐标

subproblemSize = 10;
numObjectives = M; %目标函数个数
numSubproblems = ceil(NP / subproblemSize);
F0 = 1;
CR0 = 0.5;


%% 主循环
for times = 1:num_experiments
    rng(times);
    %% 输入数据
    % 输入买家和卖家的数量
    % m = 60+5*(times-1); % 买家数量
    m = 60;
    gbest = zeros(1,m*2+2);
    n = 16; % 卖家数量
    swt = 0; % 1：添加紧急程度；0：取消紧急程度
    rep = 0; % 1:存在信誉变化；0：不存在信誉变化
    % 调用 generate_data 函数生成需求和供给数据
    [com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep);
    dim = m * 2 ;      % 决策变量维数
    maxGenerations = maxgen_base;
    % 初始化种群
    population = initPop(NP, numObjectives, dim, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);    
    
    % 生成权重向量
    weightVectors = generateWeightVectors(numSubproblems,numObjectives);

    % 非支配排序（non-domination-sort）
    population = nonDominatedSort(population, numObjectives, dim);
    start_time = tic;
    elapsed_time = zeros(maxGenerations, 1);
    for gen = 1:maxGenerations
        tic; % 记录迭代开始时间
        times
        gen
        for i = 1:numSubproblems
            lamda = exp(1-maxGenerations/(maxGenerations+1-gen));
            F = F0*2^lamda;
            CR = CR0*(1-lamda);
            % 根据权重向量和种群，选择邻居个体
            neighbors = selectNeighbors(i, weightVectors, population, subproblemSize); 
            % 交叉和变异操作
%             offspring = geneOperator(neighbors, numObjectives, dim, F, CR, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p); 
            % 执行子问题上的进化操作（交叉、变异等）
            offspring = evolveSubproblem(gbest, i, neighbors(:,1:dim+2), weightVectors, F, CR, m, n, com, spc, COM, SPC, N, r, Ur, D, x0, rho, v, ka, epsilon, sigma, p);
            % 更新外部存档
%             population = updateExternalArchive(population, offspring);
            Nc = size(population,1);
            Nn = size(offspring,1);
            allchrom(1:Nc,:) = population;
            allchrom(Nc + 1 : Nc + Nn,1 : numObjectives+dim) = offspring;
            allchrom = nonDominatedSort(allchrom, numObjectives, dim);
            population = replace_chromosome(allchrom, numObjectives, dim, NP);
            FG1(gen,1) = -min(population(:,dim + 1));
            FG2(gen,1) = min(population(:,dim + 2));
            figure(1);
            plot(-population(:,dim + 1) ,population(:,dim + 2),'*');
            str = sprintf('多目标遗传算法帕累托求解第%d次迭代',gen);
            title(str)
            xlabel('卖价与买价差值');
            ylabel('资源耗能和');
            pause(0.05)

            % 提取当前迭代的Pareto前沿

        end
            gbest = population(:,1:dim+numObjectives);
            gbest = unique(gbest,'rows');
            % 更新权重向量
            weightVectors = updateWeightVectors(weightVectors, gbest);
        % 保存当前迭代的结果
%         saveResults(paretoFront, times, gen);
            toc
            delay(gen) = toc; % 记录迭代结束时间，并计算出延迟
            %% 计算并打印累积用时
            elapsed_time(gen) = toc(start_time);
            fprintf('累积用时：%.2f秒\n', toc(start_time));
    end
    %% 平均延迟
    delay_average_results(times) = elapsed_time(maxGenerations)/maxGenerations;
        figure(1)
        plot(-population(:,dim + 1),population(:,dim + 2),'ko')
        xlabel('市场总价值')
        ylabel('资源耗能和')
        grid on
        title('帕累托解集')
        figure(2)
        plot(FG1,'k-')
        xlabel('迭代次数')
        ylabel('市场总价值')
        grid on
        
        figure(3)
        plot(FG2,'k-')
        xlabel('迭代次数')
        ylabel('资源耗能和')
        grid on
        fprintf('帕累托解集数量为 %d\n',size(gbest,1))

        if swt == 0
            file_path = '../MODEA_results'; % 修改为你希望保存的文件夹路径
        elseif swt == 1
            file_path = '../MODEA_results_Ur'; % 修改为你希望保存的文件夹路径
        end
        %保存帕累托解集
        file_name_01 = 'Pareto_results.xlsx'; % 修改为你希望保存的文件名
        gbest(:, 2*m+1) = -gbest(:, 2*m+1);   % 21列取反
        file_restore_01 = gbest(:,[2*m+1:2*m+2]);   %保存21和22列
        
        %保存总价值
        file_name_02 = 'Revenue_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_02 = FG1;   %保存21和22列
        
        %保存总能耗
        file_name_03 = 'Consumption_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_03 = FG2;   %保存21和22列

        %保存算法平均延迟
        file_name_05 = 'delay_average_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_05 = delay_average_results;

        %保存算法总延迟
        file_name_06 = 'elapsed_time_results.xlsx'; % 修改为你希望保存的文件名
        file_restore_06 = elapsed_time;
        % 使用 xlswrite 函数保存数据到 Excel 文件中
        xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_05), file_restore_05, times, 'A1'); % 将数据从 A1 单元格开始保存
        xlswrite(fullfile(file_path, file_name_06), file_restore_06, times, 'A1'); % 将数据从 A1 单元格开始保存
end
