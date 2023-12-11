function fit = fitness(ch, m, n, com, spc, COM, SPC, N, r ,Ur, D, x0 , rho, v, ka, epsilon, sigma, p)
    x = ch(1:m);
    y = ch(m+1:m*2);
    lambda = zeros(m, n); %������Դ����
    for i = 1 : m
       lambda(i,x(i)) = 1; 
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
    
    R = zeros(1, n);
    for j = 1 : n
       R(j) = 0;
       rRSU(j)=0;
       for i = 1 : m
          R(j) = R(j) +  (N(j) * r(i,j)* Ur(i) - 0.5) * (lambda(i,j) * com(i) + u(i,j) * spc(i)); % ��j��RSU��������
          rRSU(j) = r(i,j)+ rRSU(j);
       end
    end
    
    u = zeros(1, n);
    cs = zeros(1, n);
    ss = zeros(1, n);
    % ����ʵ�ʹ�Ӧ��
    for k = 1 : n
        cs(k) = sum(com(x == k));
        ss(k) = sum(spc(y == k));
        u(k)=rRSU(k)*(COM(k)-cs(k)+ SPC(k)-ss(k))+u(k);   
    end
    U = ka *(-rho * sum(C) + v * sum(R))+ 1/ka *sum(u) ;
    
    %�ܺļ���
    E = epsilon * sum(com .* x0);
    for i = 1 : m
        h = 1e-8 + exp(2-5*log10(D(i,y(i))));
        E = E + p * x0(i) / (spc(i) * log2(1 + p * h /sigma^2));
    end
    
    %% Լ������
    gfun = 0;
    cs = zeros(1, n);  % ����ʵ�ʹ�Ӧ��
    for k = 1 : n
        cs(k) = sum(com(x == k));
        gfun = gfun + max(cs(k)-COM(k),0);
    end
    
    
    ss = zeros(1, n);  % ����ʵ�ʹ�Ӧ��
    for k = 1 : n
        ss(k) = sum(spc(y == k));
        gfun = gfun + max(ss(k)-SPC(k),0);
    end
    
    
    
    
    
    fit = [-U E] + gfun * 1e8;
end