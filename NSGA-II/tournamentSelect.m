function f = tournamentSelect(chrom, poolsize, toursize)
%%  ������ѡ����� 
[NP, D] = size(chrom);
rankIndex = D - 1;
distanceIndex = D;
for i = 1 : poolsize
    % ���ѡ��toursize������
    candidate = ones(1,toursize);
    for j = 1 : toursize
        % ���ѡ��һ������
        candidate(j) = randi(NP);
        if j > 1
            % ȷ����ѡ����ͬ�ĸ���
            while ~isempty(find(candidate(1 : j - 1) == candidate(j)))
                candidate(j) = randi(NP);
            end
        end
    end
    % ��ѡ������Ӧ��rankֵ��ӵ����ֵ.
    for j = 1 : toursize
        rank(j) = chrom(candidate(j),rankIndex);
        distance(j) = chrom(candidate(j),distanceIndex);
    end
    index1 = find(rank == min(rank));
    if length(index1) ~= 1
        index2 = find(distance(index1) == max(distance(index1)));
        if length(index2) ~= 1
            index2 = index2(1);
        end
        f(i,:) = chrom(candidate(index1(index2)),:);
    else
        f(i,:) = chrom(candidate(index1(1)),:);
    end
end
