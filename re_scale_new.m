function [data] = re_scale_new(data)
    [Nb_s Nb_b]=size(data);
    M=max(data,[],1);
    m=min(data,[],1);
    data = (data-repmat(m,Nb_s,1))./(repmat(M-m,Nb_s,1));
end
