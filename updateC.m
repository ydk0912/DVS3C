function [U] = updateC(A,m,w,num)
islocal=1;
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
        q(v,:) = ai/sum(w);
    end
    U(i,idxs0) = SloutionToP19(q,m);
    clear q;
end
end

