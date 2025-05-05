function [W]= DVS3C(Nc,lamada,alpha,img)
%% 视图构造
%原始视图
[no_lines, no_rows,no_bands] = size(img);
X1=reshape(img,no_lines*no_rows, no_bands);
[N, L] = size(X1);
%超像素视图
num_Pixel=fix(N/(Nc^2));
[data2,~]=SuperpixelView(img,num_Pixel);
X2=reshape(data2,no_lines*no_rows, no_bands);
%归一化
X1=re_scale_new(X1);
X2=re_scale_new(X2);

%% 参数初始化
max_iter = 60  ;
I = eye(L);
m=2;
umax=1e2;
p=1.5;
u1 = 0.01;
u2 = 0.01;
u={u1,u2};
X={X1,X2};
%系数矩阵
A1= zeros(L);A2 = zeros(L);A={A1,A2};
J1= zeros(L);J2 = zeros(L);J={J1,J2};
M1=zeros(L);M2=zeros(L);M={M1,M2};
%权重系数
wv = ones(m, 1) / m;
%共识矩阵
C = zeros(L);
%结构相似矩阵
Sv1=zeros(L);Sv2=zeros(L);S={Sv1,Sv2};
%初始化S
[Sv1, ~] = InitializeSIGs(X1);
[Sv2, ~] = InitializeSIGs(X2);
%初始化Ls
Ls1 = diag(sum(Sv1))-Sv1;Ls2 = diag(sum(Sv2))-Sv2;Ls={Ls1,Ls2};
obj1=[];obj2=[];objv={obj1,obj2};

%% 迭代更新
for iter = 1:max_iter
    for v=1:m
        Av=A{v};
        uv=u{v};
        Mv=M{v};
        %更新J
        [ U1, S1, V1] = svd(Av+(1/uv).*Mv,'econ');
        S_thresholded =max(S1-(lamada/uv).*I,0);
        Jv=U1* S_thresholded * (V1)';
        J{v}=Jv;
        clear U1 S1 V1 Jv ;
    end

    %更新 S
    for v = 1:m
        Av=A{v};
        distU = L2_distance_1(Av,Av);
        Sv = zeros(L);
        hv = -(distU)/(2*(alpha));
        hv = (hv + hv')/2;
        for b = 1:L
            Sv(b, :) = EProjSimplex_new(hv(b, :));
        end
        Sv = (Sv+Sv')/2;
        S{v} = Sv;
        %更新L

        Ls{v} = diag(sum(S{v}))- S{v};
        clear distU;
    end

    %更新 A
    for v=1:m
        Xv = X{v};
        uv=u{v};
        Mv=M{v};
        Jv=J{v};
        A_new=pinv((Xv')*Xv+uv.*I+2.*Ls{v}+2*wv(v).*I)*((Xv')*Xv-Mv+uv.*Jv+2*wv(v).*C);
        A_new(A_new<0)=0;
        A_new=A_new-diag(diag(A_new));
        A{v}=A_new;
        %%更新M
        M{v}=Mv+uv.*(A{v}-Jv);
        %% 更新u
        u{v}=min(p*uv,umax);
    end

    %更新w{v}
    for v = 1:m
        US = C - A{v};
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end
        wv(v) = 0.5/sqrt(distUS);
    end

    %更新C
    C=updateC(A,m,wv,L);
    C= (C + C')/2;
    for v=1:m
        Xv=X{v};
        Av=A{v};
        [~, SA, ~] = svd(Av,'econ');
        temp1{v}(iter)=norm(Xv-Xv*Av,'fro')^2;
        temp2{v}(iter) = lamada*sum(sum(SA));
        temp3{v}(iter)=alpha*norm(S{v},'fro')^2+sum(sum(L2_distance_1(A{v},A{v}).* S{v}));
        temp4{v}(iter)=wv(v)*norm(C-Av,'fro')^2;
        objv{v}(iter)=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          temp1{v}(iter)+temp2{v}(iter)+temp3{v}(iter)+temp4{v}(iter);
    end
    obj1=objv{1};
    obj2=objv{2};
    obj=obj1+obj2;

    %停止条件
    if  (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6 )
        break;
    end
    if iter == 1 || mod(iter, 10) == 0
        fprintf('%f iteration compelte, obj=%f\n',iter, objv{1}(iter)+objv{2}(iter));
    end
end
W=0.5*(abs(C)+abs(C'));
end
