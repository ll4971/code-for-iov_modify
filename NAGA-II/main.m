clc
clear
close all
%% ��������
com = [30	33	36	46	49	50	54	57	59	59];   % ������������Դ
spc = [20	20	22	27	30	30	32	33	34	46];   % �������Ƶ����Դ
COM = [80	50	100	70	150	200	90];   % ���ҹ���������Դ
SPC = [40	60	80	60	100	140	100];  % ���ҹ���������Դ
N = [160	150	180	165	190	200	170];  % ����
r = [0.7	0.2	0.3	0.8	0.1	0.8	0.7	0.3	1	0.7
    0.6	0.2	0.3	0.7	0.2	0.9	0.8	0.4	0.8	0.6
    0.6	0.3	0.5	0.6	0.1	0.8	0.8	0.3	0.9	0.6
    0.6	0.4	0.4	0.7	0.1	0.7	0.7	0.5	0.9	0.6
    0.7	0.4	0.2	0.8	0.1	0.8	0.9	0.4	0.9	0.5
    0.5	0.4	0.2	0.9	0.2	0.7	0.8	0.5	0.9	0.7
    0.6	0.2	0.3	0.6	0.3	0.9	0.7	0.3	0.8	0.6
    ] * 10; %����ֵ
r = r'; %�д���V���д���RSU��[1,2]����V1��RSU2
D = [15	18	25	27	30	17	10	18	7	5
    20	13	17	19	5	7	17	30	13	12
    12	30	12	24	17	19	23	2	7	17
    8	20	20	5	14	15	29	10	18	7
    23	8	15	21	27	29	16	14	14	8
    11	25	18	16	28	9	23	23	11	17
    12	12	16	23	28	14	17	14	77	10
    ];  %����
D = D';
x0 = [150  160 185 220	222	230	250	275	290	300];
rho = 0.5; % �������Ŀ�꺯����һ���ռ��
ka = 300;   % 
v = 0.5;    % �������Ŀ�꺯����������ռ��
epsilon = 1;
sigma = 1;
p = 0.5;
m = length(com);  % �������
n = length(COM);  % ��������
%% �㷨����
NP = 80;          % ��Ⱥ����
maxgen = 600;     % ��������
Pc = 0.8;
Pm = 0.2;
M = 2;            % Ŀ�꺯������
dim = m * 2 ;      % ���߱���ά��
%% ��ʼ����Ⱥ
chrom = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,D, x0 ,rho, v,ka, epsilon, sigma, p);

%% ��֧������non-domination-sort��
chrom = nonDominatedSort(chrom, M, dim );

%% ��������
figure(1);
for gen = 1 : maxgen
    gen
    % ѡ�񸸴����ڷ�ֳ�����ԭʼNSGA-II���û���ӵ���ȱȽ����ӵĶ����ƽ���������ѡ���
    pool = round(NP/2);
    tour = 2;
    % ѡ�����    
    parentchrom = tournamentSelect(chrom, pool, tour);   
    % ����ͱ������
    newchrom = geneOperator(parentchrom, M, dim, Pm, Pc, m, n, com, spc, COM, SPC, N, r ,D, x0 , rho, v,ka, epsilon, sigma, p);   
    % �Ӵ�����͸��������ں�
    Nc = size(chrom,1);
    Nn = size(newchrom,1);
    allchrom(1:Nc,:) = chrom;
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