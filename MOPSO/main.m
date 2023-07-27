clc
clear
close all

num_experiments = 10;
distanceGenerate(); % ��þ���
for times = 1:num_experiments
    rng(times);
    %% ��������
    com = [30	33	36	46	49	50	54	57	59	59];   % ������������Դ
    x0 = [150 168 182 230 245 248 266 279 302 311]; % �������������
    spc = [20	20	22	27	30	30	32	33	34	46];   % �������Ƶ����Դ
    COM = [80	50	100	70	150	200	90];   % ���ҹ���������Դ
    SPC = [40	60	80	60	100	140	100];  % ���ҹ���������Դ
    N = [160	150	180	165	190	200	170];  % ����
    
    swt = 0; % 1����ӽ����̶ȣ�0��ȡ�������̶�
    if swt == 1
        Ur = [0.1 0.6 0.1 0.1 0.6 0.1 0.1 1 0.6 0.1]*5; %�����̶�
    elseif swt ==0
        Ur = ones(1,10);
    end
    
    r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	1	0.7
        0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.8	0.6
        0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.9	0.6
        0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.9	0.6
        0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.9	0.5
        0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.9	0.7
        0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.8	0.6
        ] * 10;
    r = r';

    rho = 0.5; % �������Ŀ�꺯����һ���ռ��
    ka = 10;   %
    v = 0.5;    % �������Ŀ�꺯����������ռ��
    epsilon = 1;
    sigma = 0.8;
    p = 1000; % ���书��
    m = length(com);  % �������
    n = length(COM);  % ��������

    dim = m * 2;
    % MOPSO
    %% ����Ⱥ����
    NP = 80;                   % ��Ⱥ����
    M = 2;
    maxgen = 600;               % ����������
    w = 0.4;                    % ����ϵ��
    c1 = 2;
    c2 = 2;
    ngrid = 20;
    Xmin = ones(1, dim);
    Xmax = n * ones(1, dim);
    Vmax = (Xmax - Xmin) / 2;
    Vmin = - Vmax;
    %% ��ʼ��
    X = zeros(NP, dim);
    V = zeros(NP,dim);
    fx = zeros(NP, M);
    for i = 1:NP
        X(i,:) = initpop(m, n, com, spc, COM, SPC);
        V(i,:) = Vmin + (Vmax - Vmin) .* rand(1,dim);
        fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v, ka, epsilon, sigma, p);
    end
    
    pbest = X;
    fpbest = fx;
    Idx = getNondominated(fx);
    rep.X  = X(Idx,:);         % ���ӽ⼯
    rep.fx = fx(Idx,:);
    rep = updateGrid(rep,ngrid);
    figure(1)
    for gen = 1:maxgen
        times
        gen
        % ѡ��leader
        h = selectLeader(rep);
        gbest = rep.X(h,:);
        for i = 1 : NP
            % �����ٶ���λ�ø���
            V(i,:) = w .* V(i,:) + c1 * rand(1,dim) .* (pbest(i,:) - X(i,:)) + c2 * rand(1,dim) .* (gbest - X(i,:));
            % ��֤�����ٶ�λ�ڽ���
            index = (V(i,:) > Vmax);
            V(i,index) = Vmax(index);
            index = (V(i,:) < Vmin);
            V(i,index) = Vmin(index);
            
            % ��Ⱥ����
            X(i,:) = X(i,:)+V(i,:);
            % ��֤����λ���ڽ���
            index = (X(i,:) > Xmax);
            X(i,index) = Xmax(index);
            index = (X(i,:) < Xmin);
            X(i,index) = Xmin(index);
            
            % ������Ӧ��
            fx(i,:) = fitness(X(i,:), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);
        end
        % ����repository
        rep = updateRepository(rep,X,fx,ngrid);
        if(size(rep.X,1) > NP)
            rep = deleteFromRepository(rep,size(rep.X,1)-NP,ngrid);
        end
        
        % ���¸�������pbest
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
        % ��ͼ��ʾ
        plot(-rep.fx(:,1),rep.fx(:,2),'ok'); hold on;
        hold off
        grid on; xlabel('f1'); ylabel('f2');
        pause(0.01)
        drawnow;
        axis square;
    end
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
    fprintf('�����н⼯����Ϊ %d\n',size(gbest,1))
    %%������excel
    file_path = '../MOPSO_results'; % �޸�Ϊ��ϣ��������ļ���·��

    %���������н⼯
    file_name_01 = 'Pareto_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_01 = fgbest;   
    
    %�����ܼ�ֵ
    file_name_02 = 'Revenue_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_02 = FG1;   
    %�������ܺ�
    file_name_03 = 'Consumption_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_03 = FG2;   

    % ʹ�� xlswrite �����������ݵ� Excel �ļ���
    xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
end