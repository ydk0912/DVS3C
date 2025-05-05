function v = ToVector(im)
% takes MxNx3 picture and returns (MN)x3 vector
% im m*n*b
sz = size(im);   % 返回sz = [m,n,b] 145,145,200
v = reshape(im, [prod(sz(1:2)) sz(3)]); %都是将A 的行列排列成m行n列。另外 reshape是 按照列取数据的
                                        %prod() 列的元素连乘乘积 prod（sz（1:2）） =145*145
                                        %v = [145*145,200]
