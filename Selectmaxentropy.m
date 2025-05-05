function [Y] = Selectmaxentropy(labels, img)
% %     %用类中熵最大的波段作为该类被选择的波段
    num_clusters = max(labels); % 类簇的数量
    Y=[];
    for i = 1:num_clusters
        cluster_points = find(labels == i); % 获取第i个类簇的样本
        p_num=length(cluster_points);
        entropy_d=[];
        for j=1:p_num
            entropy_d(j)=entropy(img(:,:,cluster_points(j)));
        end
        [~, center_index] = max(entropy_d); % 
        Y =[Y,cluster_points(center_index)]; 
    end