function v = ToVector(im)
% takes MxNx3 picture and returns (MN)x3 vector
% im m*n*b
sz = size(im);   % ����sz = [m,n,b] 145,145,200
v = reshape(im, [prod(sz(1:2)) sz(3)]); %���ǽ�A ���������г�m��n�С����� reshape�� ������ȡ���ݵ�
                                        %prod() �е�Ԫ�����˳˻� prod��sz��1:2���� =145*145
                                        %v = [145*145,200]
