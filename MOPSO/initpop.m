function ch = initpop(m, n, com, spc, COM, SPC)

% ��ʼ����Ⱥ
while 1
    x = randi(n, 1, m);
    cs = zeros(1, n);  % ����ʵ�ʹ�Ӧ��
    for k = 1 : n
        cs(k) = sum(com(x == k));
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
ch = [x y];

