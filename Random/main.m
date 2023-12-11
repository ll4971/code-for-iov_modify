%% ���ö��ʵ��
clc
clear
close all
num_experiments = 1;
%% ��������
% ������Һ����ҵ�����
m = 50; % �������
n = 7; % ��������
swt = 1; % 1����ӽ����̶ȣ�0��ȡ�������̶�
rep = 0; % 1:���������仯��0�������������仯
% ���� generate_data ������������͹�������
[com, spc, COM, SPC, Ur, r, N, D, x0] = generate_data(m, n, swt, rep);


    rho = 0.5; % �������Ŀ�꺯����һ���ռ��
    ka = 10;   %
    v = 0.5;    % �������Ŀ�꺯����������ռ��
    epsilon = 1;
    sigma = 6;
    p = 1000; % ���书��
    %% �㷨����
    NP = 80;          % ��Ⱥ����
    maxgen = 60;     % ��������
    Pc = 0.8;
    Pm = 0.2;
    M = 2;            % Ŀ�꺯������
    dim = m * 2 ;      % ���߱���ά��
%% ��ѭ��
for times = 1:num_experiments
    rng(8);

    %% ��¼�ӳٵ�����
    delay = zeros(maxgen, 1);
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    figure(1);
    for gen = 1 : maxgen
        tic; % ��¼������ʼʱ��
        times
        gen
        %% ��ʼ����Ⱥ
        chrom = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 ,rho, v,ka, epsilon, sigma, p);
    
        %% ��֧������non-domination-sort��
        chrom = nonDominatedSort(chrom, M, dim );
        FG1(gen,1) = -min(chrom(:,dim + 1));
        FG2(gen,1) = min(chrom(:,dim + 2));
        plot(-chrom(:,dim + 1) ,chrom(:,dim + 2),'*');
        str = sprintf('��Ŀ���Ŵ��㷨����������%d�ε���',gen);
        title(str)
        xlabel('��������۲�ֵ');
        ylabel('��Դ���ܺ�');
        pause(0.05)
        hold off
        toc
        delay(gen) = toc; % ��¼��������ʱ�䣬��������ӳ�
        %% ���㲢��ӡ�ۻ���ʱ
        elapsed_time(gen) = toc(start_time);
        fprintf('�ۻ���ʱ��%.2f��\n', toc(start_time)); 
    end
    %% ������
    clc
    gbest = chrom(:,1:dim+M);
    gbest = unique(gbest,'rows');
    
    figure(1)
    plot(-chrom(:,dim + 1),chrom(:,dim + 2),'ko')
    xlabel('�г��ܼ�ֵ')
    ylabel('��Դ���ܺ�')
    grid on
    title('�����н⼯')
    figure(2)
    plot(FG1,'k-')
    xlabel('��������')
    ylabel('�г��ܼ�ֵ')
    grid on
    figure(3)
    plot(FG2,'k-')
    xlabel('��������')
    ylabel('��Դ���ܺ�')
    grid on
    fprintf('�����н⼯����Ϊ %d\n',size(gbest,1))
    figure(4) 
    plot(elapsed_time) 
    xlabel('��������') 
    ylabel('�ۻ���ʱ') 
    title('�ۻ���ʱ����������ı仯') 
    
    %%������excel
    if swt == 0
        file_path = '../NAGA-II_results'; % �޸�Ϊ��ϣ��������ļ���·��
    elseif swt == 1
        file_path = '../NAGA-II_results_Ur'; % �޸�Ϊ��ϣ��������ļ���·��
    end
    %���������н⼯
    file_name_01 = 'Pareto_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    gbest(:, 21) = -gbest(:, 21);   % 21��ȡ��
    file_restore_01 = gbest(:,[21:22]);   %����21��22��
    
    %�����ܼ�ֵ
    file_name_02 = 'Revenue_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_02 = FG1;   %����21��22��
    
    %�������ܺ�
    file_name_03 = 'Consumption_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_03 = FG2;   %����21��22��
    
    % ʹ�� xlswrite �����������ݵ� Excel �ļ���
%     xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
%     xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
%     xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
end

