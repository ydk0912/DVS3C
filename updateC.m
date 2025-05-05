function [U] = updateC(A,m,w,num)
islocal=1;
% lambda=1;
% dist = L2_distance_1(F', F');
U = zeros(num);
for i=1:num
    idx = zeros();
    for v = 1:m
        a0 = A{v}(i,:);
        idx = [idx,find(a0>0)];
    end
    idxs = unique(idx(2:end));
    if islocal == 1
        idxs0 = idxs;
    else
        idxs0 = 1:num;
    end
    for v = 1:m
        a1 = A{v}(i,:)*w(v);
        ai = a1(idxs0);
%         di = dist(i,idxs0);
%         mw = m*w(v);
%         lmw = lambda/mw;
%         q(v,:) = ai-0.5*lmw*di;
        q(v,:) = ai/sum(w);
    end
    U(i,idxs0) = SloutionToP19(q,m);
    clear q;
end
end

%   function C = updateC(AV, num_view,Wv, num_samp)
%         C = zeros(num_samp);
%         for i = 1:num_samp
%             zv0 = zeros(1, num_samp);
%             for v = 1:num_view
%                 temp = AV{v};
%                 zv0 = zv0 + Wv(v)*temp(i,:);
%             end
%             idxa0 = find(zv0>0);
%             zi = zv0(idxa0);
% %             ui = dist(i, idxa0);
%             cu = (zi)/sum(Wv);
%             C(i,idxa0) = EProjSimplex_new(cu);
%         end
%     end