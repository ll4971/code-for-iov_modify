% ���±����� 
function rep = updateRepository(rep, X, fx, ngrid)
% ȷ����֧���
Idx  = getNondominated(fx);
rep.X    = [rep.X; X(Idx,:)];
rep.fx = [rep.fx; fx(Idx,:)];
% �ϲ���ķ�֧���
Idx  = getNondominated(rep.fx);
rep.fx= rep.fx(Idx,:);
rep.X    = rep.X(Idx,:);
% Updating the grid
rep        = updateGrid(rep,ngrid);
end
