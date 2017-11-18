function mls
% 脚本文件：moving least squares method

% 文件名：mls

% 清除空间变量
clear;
% 以下为参数的设置
% ----------------------------- %
% 已知点的坐标矩阵（矩形）
a=-1:0.01:1;
b=-1:0.01:1;
[X,Y]=meshgrid(a,b);
% 已知点的值矩阵
f=sinh(5*X)+20*sin(2*pi*Y);
% 局部支撑域的半径尺度
rs=0.025;
% 权函数的形状参数
alpha=2.5;
% 径向基函数次数(只能填1or2)
k=1;
% ----------------------------- %
c=-1:0.05:1;
d=-1:0.05:1;
[J,K]=meshgrid(c,d);
n_s=numel(J);
Z=zeros(size(J));
for i=1:n_s
    [Z(i)]=mlss(J(i),K(i),f,X,Y,alpha,rs,k);
end
surf(J,K,Z);

end
% x,y为计算点，f,X,Y为已知点参数，alpha为权函数的形状参数，rs为局部支撑域半径，k为径向基最高次数
function [ux] = mlss(x,y,f,X,Y,alpha,rs,k)
% n_p径向基元素个数，n已知散点个数
n_p=(k+1)*(k+2)/2;
n=numel(X);
% A矩阵模型,B矩阵模型，自然对数e
A=zeros(n_p,n_p);
B=zeros(n_p,n);
e=exp(1);
% 移动最小二乘核心
for i=1:n
    [p]=Radial_basis(X(i),Y(i),k);
    % Gauss型权函数部分
    r=sqrt((X(i)-x).^2+(Y(i)-y).^2)/rs;
    if r<=1
        w=e.^(alpha.*r);
    else
        w=0;
    end
    A=A+w.*(p*p');
    bb=w.*p;
    B(:,i)=bb;
    u(i,:)=f(i);
end
D=A^-1;
ax=D*B*u;
[pp]=Radial_basis(x,y,k);
ux=pp'*ax;
end

% xi,yi为已知散点坐标，k为径向基次数
function [p]=Radial_basis(xi,yi,k)
if k==1
    p=[1;xi;yi];
elseif k==2
    p=[1;xi;yi;xi*xi;xi*yi;yi*yi];
end
end


% x,y为计算点，xi,yi为已知散点坐标，k为径向基次数
function [p]=Radial_basis_kernel(xi,yi,x,y,k)
if k==1
    p=[1,xi-x,yi-y];
elseif k==2
    p=[1,xi-x,yi-y,(xi-x)*(xi-x),(xi-x)*(yi-y),(yi-y)*(yi-y)];
end
end