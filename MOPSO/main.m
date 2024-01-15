%% ���ö��ʵ��
clc 
clear
close all
num_experiments = 1;
%% �㷨����
NP = 200;          % ��Ⱥ����
maxgen_base = 1500;     % ��������
Pc = 0.8;
Pm = 0.2;
M = 2;            % Ŀ�꺯������
rho = 0.5; % �������Ŀ�꺯����һ���ռ��
ka = 10;  
v = 0.5;    % �������Ŀ�꺯����������ռ��
epsilon = 1;
sigma = 174;
omiga = 40; %����
p = 100; % ���书��
T_delay = zeros(num_experiments,1); % ���ӳ�ʱ��
delay_average_results = zeros(num_experiments,1); % ƽ���ӳ�ʱ��
x = zeros(num_experiments,1); % ������������

%% ��ѭ��
for times = 1:num_experiments
    rng(times);
    %% ��������
    % ������Һ����ҵ�����
    % m = 60+5*(times-1); % �������
    m = 60;
    n = 16; % ��������
    swt = 0; % 1����ӽ����̶ȣ�0��ȡ�������̶�
    rep1 = 0; % 1:���������仯��0�������������仯
    % ���� generate_data ������������͹�������
    [com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep1);
    dim = m * 2 ;      % ���߱���ά��
    maxgen = maxgen_base + 10*m;
    % MOPSO
    %% ����Ⱥ����
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
    rep.X  = X(Idx,:);         % ���ӽ⼯
    rep.fx = fx(Idx,:);
    rep = updateGrid(rep,ngrid);
    figure(1)
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    for gen = 1:maxgen
        tic; % ��¼������ʼʱ��
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
        toc
        delay(gen) = toc; % ��¼��������ʱ�䣬��������ӳ�
        %% ���㲢��ӡ�ۻ���ʱ
        elapsed_time(gen) = toc(start_time);
        fprintf('�ۻ���ʱ��%.2f��\n', toc(start_time));
    end
    %% ƽ���ӳ�
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

    %�����㷨ƽ���ӳ�
    file_name_05 = 'delay_average_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_05 = delay_average_results;

    %�����㷨���ӳ�
    file_name_06 = 'elapsed_time_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_06 = elapsed_time;

    % ʹ�� xlswrite �����������ݵ� Excel �ļ���
    xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_05), file_restore_05, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
    xlswrite(fullfile(file_path, file_name_06), file_restore_06, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
end