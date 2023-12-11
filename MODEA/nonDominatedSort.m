function f = nonDominatedSort(x, M, dim)
% �ú����ǶԵ�ǰ��Ⱥ���з�֧�����򡣵�һ��(Front)�ĸ����Ϊrank 1���ڶ����
% ���Ӧrank 2���������ơ�ָ��rank�󣬼���ÿһ���ӵ����. 
N = size(x,1);
front = 1;    % ��ǰ�����Ϊ��һ��.
F(front).f = [];
individual = [];
%% ��֧������
for i = 1 : N
    individual(i).np = 0;   % ֧��ø��� i �ĸ�������
    individual(i).Sp = [];  % ���ø���֧��ĸ��弯��
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M      % ��ÿһ��Ŀ�꺯��ֵ���бȽ�
            if (x(i,dim + k) < x(j,dim + k))
                dom_less = dom_less + 1;
            elseif (x(i,dim + k) == x(j,dim + k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M    % �����������Ӧ����ȥ���������˼������i~=j��Ӧ����������
            individual(i).np = individual(i).np + 1;   % q ֧�� p
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).Sp = [individual(i).Sp j];   % p ֧�� q
        end
    end   
    if individual(i).np == 0
        x(i,M + dim + 1) = 1;          % �����д��Ⱦɫ������һλ
        F(front).f = [F(front).f i]; % ��һ����弯��
    end
end

% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)    % ��front�������и����������ѭ��
       if ~isempty(individual(F(front).f(i)).Sp)   % front���е�i�������Sp��Ϊ��
        	for j = 1 : length(individual(F(front).f(i)).Sp) % ��front����i�������Sp�����и����������ѭ��
            	individual(individual(F(front).f(i)).Sp(j)).np = individual(individual(F(front).f(i)).Sp(j)).np - 1;
        	   	if individual(individual(F(front).f(i)).Sp(j)).np == 0
               		x(individual(F(front).f(i)).Sp(j),M + dim + 1) = front + 1;
                    Q = [Q individual(F(front).f(i)).Sp(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + dim + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;

%% ����ӵ����crowding distance
for front = 1 : (length(F) - 1)
%    objective = [];
    distance = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;
    sorted_based_on_objective = [];
    for i = 1 : M
        [sorted_based_on_objective, index_of_objectives] = sort(y(:,dim + i));
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives)
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
        end
        f_max =  sorted_based_on_objective(length(index_of_objectives), dim + i);
        f_min = sorted_based_on_objective(1, dim + i);
        y(index_of_objectives(length(index_of_objectives)),M + dim + 1 + i) = Inf;
        y(index_of_objectives(1),M + dim + 1 + i) = Inf;
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1,dim + i);
            previous_obj  = sorted_based_on_objective(j - 1,dim + i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + dim + 1 + i) = Inf;
            else
                y(index_of_objectives(j),M + dim + 1 + i) =  (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        distance(:,1) = distance(:,1) + y(:,M + dim + 1 + i);
    end
    y(:,M + dim + 2) = distance;
    y = y(:,1 : M + dim + 2);
    z(previous_index:current_index,:) = y;
end
f = z();
