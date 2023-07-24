% �ӱ�������ɾ������ 
function rep = deleteFromRepository(rep, Ne, ngrid)
[n,M] = size(rep.fx);
fx = rep.fx;
CD = zeros(n, M);  % ÿ���������Ŀ�꺯���µ�ӵ����
for m = 1 : M
    [~, objfunsortidx] = sort(fx(:, m));
    fxs = fx(objfunsortidx,:);
    objfunsort_min = fxs(1, m);
    bojfunsort_max = fxs(n, m);
    CD(objfunsortidx(1), m) = inf;
    CD(objfunsortidx(n), m) = inf;
    
    for j = 2 : n - 1
        objfun_next  = fxs(j + 1, m);
        objfun_last  = fxs(j - 1, m);
        if (bojfunsort_max - objfunsort_min == 0)
            CD(objfunsortidx(j), m) = Inf;
        else
            CD(objfunsortidx(j), m) =  (objfun_next - objfun_last)/(bojfunsort_max - objfunsort_min);
        end
    end
end
CD = sum(CD,2);
% ɾ�����������
[~,del_idx] = sort(CD,'ascend');
del_idx = del_idx(1:Ne);
rep.X(del_idx,:) = [];
rep.fx(del_idx,:) = [];
rep = updateGrid(rep,ngrid);
end
