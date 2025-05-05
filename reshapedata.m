function [data]=reshapedata(data)
[Nb_s Nb_b]=size(data);

    AA = double(data);
    minv = min(AA(:));
    maxv = max(AA(:));
    data= double(AA - minv) / double(maxv - minv);

