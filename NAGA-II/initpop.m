function chrom = initpop(NP, M, dim, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v,ka, epsilon, sigma, p)
% K ������Ԫ���ܸ��������߱�����Ŀ�꺯��ֵ���һ������
K = M + dim;
%% ��ʼ��ÿ��Ⱦɫ��
chrom = zeros(NP,K);
for i = 1 : NP 
    % ��ʼ����Ⱥ
    while 1
        
        x = randi(n, 1, m); % ����m��1~n���������������Ҷ�Ӧ������
        cs = zeros(1, n);  % ����ʵ�ʹ�Ӧ��
        for k = 1 : n
            cs(k) = sum(com(x == k)); % ����k�Ĺ�Ӧ��
        end
        if all(cs <= COM)
            break
        end
    end
    while 1
        y = randi(n, 1, m);
        ss = zeros(1, n);  % ����ʵ�ʹ�Ӧ��
        for k = 1 : n
            ss(k) = sum(spc(y == k));
        end
        if all(ss <= SPC)
            break
        end
    end
    chrom(i,1:dim) = [x y];
end
for i = 1 : NP
    chrom(i,dim + 1: K) = fitness(chrom(i,:), m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v, ka, epsilon, sigma, p);
end
end
