function mls
% �ű��ļ���moving least squares method

% �ļ�����mls

% ����ռ����
clear;
% ����Ϊ����������
% ----------------------------- %
% ��֪���������󣨾��Σ�
a=-1:0.01:1;
b=-1:0.01:1;
[X,Y]=meshgrid(a,b);
% ��֪���ֵ����
f=sinh(5*X)+20*sin(2*pi*Y);
% �ֲ�֧����İ뾶�߶�
rs=0.025;
% Ȩ��������״����
alpha=2.5;
% �������������(ֻ����1or2)
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
% x,yΪ����㣬f,X,YΪ��֪�������alphaΪȨ��������״������rsΪ�ֲ�֧����뾶��kΪ�������ߴ���
function [ux] = mlss(x,y,f,X,Y,alpha,rs,k)
% n_p�����Ԫ�ظ�����n��֪ɢ�����
n_p=(k+1)*(k+2)/2;
n=numel(X);
% A����ģ��,B����ģ�ͣ���Ȼ����e
A=zeros(n_p,n_p);
B=zeros(n_p,n);
e=exp(1);
% �ƶ���С���˺���
for i=1:n
    [p]=Radial_basis(X(i),Y(i),k);
    % Gauss��Ȩ��������
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
    p=[1,xi-x,yi-y];
elseif k==2
    p=[1,xi-x,yi-y,(xi-x)*(xi-x),(xi-x)*(yi-y),(yi-y)*(yi-y)];
end
end