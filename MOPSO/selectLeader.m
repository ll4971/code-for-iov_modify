function selected = selectLeader(rep) 
% ���̶�
Pcum = cumsum(rep.quality(:,2));     % �����ۻ�����
sel_hyp = rep.quality(find(rand(1,1)*max(Pcum)<=Pcum,1,'first'),1); % ѡ��������
% �ڱ�ѡ�еĳ������������ѡ��һ��������ΪLeader
idx = 1:1:length(rep.grid_idx);
selected = idx(rep.grid_idx==sel_hyp);
selected = selected(randi(length(selected)));
end

