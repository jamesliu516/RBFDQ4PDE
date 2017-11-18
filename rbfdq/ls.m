function ls
% �ű��ļ���least squares method
% �ļ�����ls

% ����ռ����
clear;
% ����Ϊ����������
% ----------------------------- %
% ��֪���������󣨾��Σ�
a=-1:0.1:1;
b=-1:0.1:1;
[X,Y]=meshgrid(a,b);
% ��֪���ֵ����
f=sinh(5*X)+20*sin(2*pi*Y);
% �������������(ֻ����1or2)
k=2;
% ----------------------------- %
% c=-1:0.1:1;
% d=-1:0.1:1;
% [J,K]=meshgrid(c,d);
% n_s=numel(J);
% Z=zeros(size(J));
% for i=1:n_s
%     [Z(i)]=mlss(J(i),K(i),f,X,Y,k);
% end
% surf(J,K,Z);
[aaa]=mlss(-1,-1,f,X,Y,k);
%-------------------------------- %

end
% x,yΪ����㣬f,X,YΪ��֪�������alphaΪȨ��������״������rsΪ�ֲ�֧����뾶��kΪ�������ߴ���
function [ux] = mlss(x,y,f,X,Y,k)
% n_p�����Ԫ�ظ�����n��֪ɢ�����
n_p=(k+1)*(k+2)/2;
n=numel(X);
% A����ģ��,B����ģ�ͣ���Ȼ����e
A=zeros(n_p,n_p);
B=zeros(n_p,n);
e=exp(1);
% ��С���˺���
for i=1:n
    [p]=Radial_basis(X(i),Y(i),k);
    A=A+(p*p');
    bb=p;
    B(:,i)=bb;
    u(i,:)=f(i);
end 
D=A^-1;
ax=D*B*u;
[pp]=Radial_basis(x,y,k);
ux=pp'*ax;
end

% xi,yiΪ��֪ɢ�����꣬kΪ���������
function [p]=Radial_basis(xi,yi,k)
if k==1
    p=[1;xi;yi];
elseif k==2
    p=[1;xi;yi;xi*xi;xi*yi;yi*yi];
end   
end


% x,yΪ����㣬xi,yiΪ��֪ɢ�����꣬kΪ���������
function [p]=Radial_basis_kernel(xi,yi,x,y,k)
if k==1
    p=[1;xi-x;yi-y];
elseif k==2
    p=[1;xi-x;yi-y;(xi-x)*(xi-x);(xi-x)*(yi-y);(yi-y)*(yi-y)];
end   
end