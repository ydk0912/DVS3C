function [X2,sup_img] = SuperpixelView(data,Pix_num)
    img=data;
    [w,h,b]=size(data);
    labels=cubseg(data,Pix_num);
    labels=labels+1;
    [data,sup_img]=generatemean(img,labels);
    data=weighted_mean_feature(labels,img,sup_img);
    data=reshape(data,w,h,b);
    X2=data;
    sup_img=sup_img';
end