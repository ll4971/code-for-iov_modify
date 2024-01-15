%% 设置多次实验
clc 
clear
close all
num_experiments = 1;
%% 算法数据
NP = 200;          % 种群数量
maxgen_base = 1500;     % 迭代次数
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

%% 主循环
for times = 1:num_experiments
    rng(times);
    %% 输入数据
    % 输入买家和卖家的数量
    % m = 60+5*(times-1); % 买家数量
    m = 60;
    n = 16; % 卖家数量
    swt = 0; % 1：添加紧急程度；0：取消紧急程度
    rep1 = 0; % 1:存在信誉变化；0：不存在信誉变化
    % 调用 generate_data 函数生成需求和供给数据
    [com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep1);
    dim = m * 2 ;      % 决策变量维数
    maxgen = maxgen_base + 10*m;
    % MOPSO
    %% 粒子群参数
    w = 0.4;                    % 惯性系数
    c1 = 2;
    c2 = 2;
    ngrid = 20;
    Xmin = ones(1, dim);
    Xmax = n * ones(1, dim);
    Vmax = (Xmax - Xmin) / 2;
    Vmin = - Vmax;
    %% 初始化
    X = zeros(NP, dim);
    V = zeros(NP, dim);
    fx = zeros(NP, M);
    for i = 1:NP
        X(i,:) = initpop(m, n, com, spc, COM, SPC);
        V(i,:) = Vmin + (Vmax - Vmin) .* rand(1,dim);
        fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v, ka, epsilon, sigma, p);
    end
    
    pbest = X;
    fpbest = fx;
    Idx = getNondominated(fx);
    rep.X  = X(Idx,:);         % 非劣解集
    rep.fx = fx(Idx,:);
    rep = updateGrid(rep,ngrid);
    figure(1)
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    for gen = 1:maxgen
        tic; % 记录迭代开始时间
        times
        gen
        % 选择leader
        h = selectLeader(rep);
        gbest = rep.X(h,:);
        for i = 1 : NP
            % 粒子速度与位置更新
            V(i,:) = w .* V(i,:) + c1 * rand(1,dim) .* (pbest(i,:) - X(i,:)) + c2 * rand(1,dim) .* (gbest - X(i,:));
            % 保证粒子速度位于界内
            index = (V(i,:) > Vmax);
            V(i,index) = Vmax(index);
            index = (V(i,:) < Vmin);
            V(i,index) = Vmin(index);
            
            % 种群更新
            X(i,:) = X(i,:)+V(i,:);
            % 保证粒子位置在界内
            index = (X(i,:) > Xmax);
            X(i,index) = Xmax(index);
            index = (X(i,:) < Xmin);
            X(i,index) = Xmin(index);
            
            % 评价适应度
            fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);
        end
        % 更新repository
        rep = updateRepository(rep,X,fx,ngrid);
        if(size(rep.X,1) > NP)
            rep = deleteFromRepository(rep,size(rep.X,1)-NP,ngrid);
        end
        
        % 更新个体最优pbest
        if dominates(fx(i,:), fpbest(i,:))
            fpbest(i,:) = fx(i,:);
            pbest(i,:) = X(i,:);
        elseif dominates(fpbest(i,:), fx(i,:))
            
        else
            if rand < 0.5
                fpbest(i,:) = fx(i,:);
                pbest(i,:) = X(i,:);
            end
        end
        FG1(gen,1) = -min(rep.fx(:,1));
        FG2(gen,1) = min(rep.fx(:,2));
        % 绘图显示
        plot(-rep.fx(:,1),rep.fx(:,2),'ok'); hold on;
        hold off
        grid on; xlabel('f1'); ylabel('f2');
        pause(0.01)
        drawnow;
        axis square;
        toc
        delay(gen) = toc; % 记录迭代结束时间，并计算出延迟
        %% 计算并打印累积用时
        elapsed_time(gen) = toc(start_time);
        fprintf('累积用时：%.2f秒\n', toc(start_time));
    end
    %% 平均延迟
    delay_average_results(times) = elapsed_time(maxgen)/maxgen;
    gbest = rep.X;
    fgbest = rep.fx;
    [fgbest,IA] = unique(fgbest,'rows');
    gbest = round(gbest(IA,:));
    fgbest(:,1) = -fgbest(:,1);
    figure(1)
    plot(fgbest(:,1),fgbest(:,2),'ko')
    xlabel('Overall Revenue')
    ylabel('Resource Consumption')
    grid on
    title('Pareto solution set')
    figure(2)
    plot(FG1,'k-')
    xlabel('Iterations Times')
    ylabel('Overall Revenue')
    grid on
    
    figure(3)
    plot(FG2,'k-')
    xlabel('Iterations Times')
    ylabel('Resource consumption')
    grid on
    fprintf('帕累托解集数量为 %d\n',size(gbest,1))
    %%保存至excel
    file_path = '../MOPSO_results'; % 修改为你希望保存的文件夹路径

    %保存帕累托解集
    file_name_01 = 'Pareto_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_01 = fgbest;   
    
    %保存总价值
    file_name_02 = 'Revenue_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_02 = FG1;   
    %保存总能耗
    file_name_03 = 'Consumption_results.xlsx'; % 修改为你希望保存的文件名
    file_restore_03 = FG2;

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