clc;clear;
data = load("Indian_pines.mat").indian_pines;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
img=reshapedata(data);
Nc=6;
lamada=1;
alpha=1;
W=DVS3C(Nc,lamada,alpha,img);
clu_num=5:5:50;
for clu_i=1:length(clu_num)
    clu=clu_num(clu_i);
    labels=clu_ncut(W,clu);
    BS=Selectmaxentropy(labels,img);
end









































