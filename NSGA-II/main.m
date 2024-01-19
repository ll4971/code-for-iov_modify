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
omiga = 40; % ����
p = 100; % ���书��
T_delay = zeros(num_experiments,1); % ���ӳ�ʱ��
delay_average_results = zeros(num_experiments,1); % ƽ���ӳ�ʱ��
x = zeros(num_experiments,1); % ������������
%% ��ѭ��
for times = 1:num_experiments
    rng(10);
    %% ��������
    % ������Һ����ҵ�����
    m = 60+5*(times-1); % �������
    %m = 60;
    n = 16; % ��������
    swt = 0; % 1����ӽ����̶ȣ�0��ȡ�������̶�
    rep = 0; % 1:���������仯��0�������������仯
    % ���� generate_data ������������͹�������
    [com, spc, COM, SPC, Ur, r, N, D, x0, Ur_order] = generate_data(m, n, swt, rep);
    dim = m * 2 ;      % ���߱���ά��
    maxgen = maxgen_base + 10*m;
    %% ��ʼ����Ⱥ
    [chrom, com, spc] = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 ,rho, v,ka, epsilon, sigma, p);
    
    %% ��֧������non-domination-sort��
    chrom = nonDominatedSort(chrom, M, dim );
    %% ��¼�ӳٵ�����
    delay = zeros(maxgen, 1);
    start_time = tic;
    elapsed_time = zeros(maxgen, 1);
    %% ��������
    figure(1);
    for gen = 1 : maxgen
        tic; % ��¼������ʼʱ��
        times
        gen
        % ѡ�񸸴����ڷ�ֳ�����ԭʼNSGA-II���û���ӵ���ȱȽ����ӵĶ����ƽ���������ѡ���
        pool = round(NP/2);
        tour = 2;
        % ѡ�����    
        parentchrom = tournamentSelect(chrom, pool, tour);   
        % ����ͱ������
        newchrom = geneOperator(parentchrom, M, dim, Pm, Pc, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p);   
        % �Ӵ�����͸��������ں�
        Nc = size(chrom,1);
        Nn = size(newchrom,1);
        allchrom = chrom;
        allchrom(Nc + 1 : Nc + Nn,1 : M+dim) = newchrom;
      
        % �ں���Ⱥ��֧������
        allchrom = nonDominatedSort(allchrom, M, dim);
        chrom = replace_chromosome(allchrom, M, dim, NP);
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

    %% �����ӳ� + �����ӳ�
    y = chrom(1, m + 1: dim);
    T_delay(times) = T_delay(times) + sum(x0 ./ com);
    for i = 1:m
        h = 128.1 + 37.5*(log10(D(i,y(i))));
        T_delay(times) = T_delay(times) + x0(i) / (0.04* spc(i) * log2(1 + p * h /sigma^2));%40�Ǵ���
    end
    %% ƽ���ӳ�
    delay_average_results(times) = elapsed_time(maxgen)/maxgen;

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
    x(times) = m; 
    figure(5) 
    plot(x, T_delay)
    xlabel('������') 
    ylabel('���ӳ�') 
    title('���ӳ��泵�����ı仯')

    %% ����������Դ����
    xt = gbest(1,1:m);
    y = gbest(1,m+1:m*2);
    lambda = zeros(m, n); %������Դ����
    for i = 1 : m
       lambda(i,xt(i)) = 1; 
    end
    u = zeros(m, n); %Ƶ����Դ����
    for i = 1 : m
       u(i,y(i)) = 1; 
    end

    C = zeros(1, m);
    for i = 1 : m
       C(i) = 0;
       for j = 1 : n
          C(i) = C(i) +  N(j) / r(i,j) /Ur(i)* (lambda(i,j) * com(i) + u(i,j) * spc(i)); % ��i��V���򻨷�
       end
    end

    %% RSU������
    R = zeros(1, n);
    percent = zeros(1,n);
    R_1 = zeros(1, n);
    for j = 1 : n
        R(j) = 0;
        for i = 1 : m
            if Ur_order(i) == 5
                R_1(j) = R_1(j) + (N(j) * r(i,j)* Ur(i) - 5e-3) * (lambda(i,j) * com(i) + u(i,j) * spc(i));
            end
            R(j) = R(j) +  (N(j) * r(i,j)* Ur(i) - 5e-3) * (lambda(i,j) * com(i) + u(i,j) * spc(i)); % ��j��RSU��������
        end
        percent(j) = R_1(j)/R(j);
    end
    %% ������excel
    if swt == 0
        file_path = '../NSGA-II_results'; % �޸�Ϊ��ϣ��������ļ���·��
    elseif swt == 1
        file_path = '../NSGA-II_results_Ur'; % �޸�Ϊ��ϣ��������ļ���·��
    end
    %���������н⼯
    file_name_01 = 'Pareto_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    gbest(:, dim+1) = -gbest(:, dim+1);   % 21��ȡ��
    file_restore_01 = gbest(:,[dim+1:dim+2]);   %����21��22��
    
    %�����ܼ�ֵ
    file_name_02 = 'Revenue_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_02 = FG1;   %����21��22��
    
    %�������ܺ�
    file_name_03 = 'Consumption_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_03 = FG2;   %����21��22��

    %�������ӳ�
    file_name_04 = 'T_delay_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_04 = [x, T_delay];

    %�����㷨ƽ���ӳ�
    file_name_05 = 'delay_average_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_05 = delay_average_results;
    
    %�����㷨���ӳ�
    file_name_06 = 'elapsed_time_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_06 = elapsed_time;

    %���泵��������Դ����
    file_name_07 = 'resource_consumption_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_07 = C;
    
    %����RSU����
    file_name_08 = 'RSU_revenue_results.xlsx'; % �޸�Ϊ��ϣ��������ļ���
    file_restore_08 = [R, percent];

    % ʹ�� xlswrite �����������ݵ� Excel �ļ���
     xlswrite(fullfile(file_path, file_name_01), file_restore_01, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_02), file_restore_02, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_03), file_restore_03, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_04), file_restore_04, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_05), file_restore_05, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_06), file_restore_06, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_07), file_restore_07, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
     xlswrite(fullfile(file_path, file_name_08), file_restore_08, times, 'A1'); % �����ݴ� A1 ��Ԫ��ʼ����
end

