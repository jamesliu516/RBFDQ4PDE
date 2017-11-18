close all;
clc;clear;
q=1.15;
n=15;             %长度份数
m=16;             %角度份数
width=5;
radmin=1;
r=3;
%----------------------------------------
num=(n-r).*m/2;
x1=zeros(1,num);
y1=zeros(1,num);
z=0;
for i=1:m/2
    for j=1:n-r
        z=z+1;
        a1=(width/cos(pi.*(i-1)./(2*m))-radmin)./(1-q.^n)*(1-q);
        x1(1,z)=(a1.*(1-q.^(j+r))/(1-q)+radmin).*cos(pi.*(i-1)./(2*m));
        y1(1,z)=(a1.*(1-q.^(j+r))/(1-q)+radmin).*sin(pi.*(i-1)./(2*m));
    end
end
x2=zeros(1,n-r);
% y2=zeros(1,n);
a2=(width/cos(pi.*(m/2)./(2*m))-radmin)./(1-q.^n)*(1-q);
for k=1:n-r
    x2(1,k)=(a2.*(1-q.^(k+r))/(1-q)+radmin).*cos(pi.*(m./2)./(2*m));
%     y2(1,k)=(a2.*(1-q.^k)/(1-q)+radmin).*sin(pi.*(m./2)./(2*m));
end
x3=zeros(1,m);
y3=zeros(1,m);
ff=0;
% s3=a2.*(1-q.^(r+1))/(1-q);
% x4=((s3+radmin).*cos(pi.*(m./2)./(2*m))+x2(1,1))/2;
% y4=((s3+radmin).*sin(pi.*(m./2)./(2*m))+y2(1,1))/2;
a3=(width-radmin)./(1-q.^n)*(1-q);
sr=zeros(1,r+1);
for v=1:r+1
    sr(1,v)=a3.*(1-q.^(v-1))/(1-q);
end
for g=1:r+1
    for f=1:m+1
        ff=ff+1;
        rr=sr(1,g);
        x3(1,ff)=(rr+radmin).*cos(pi.*((f-1)./2)./m);
        y3(1,ff)=(rr+radmin).*sin(pi.*((f-1)./2)./m);
    end
end
% s4=((r-1)*(s3/(r))+radmin).*cos(pi.*((m/2)./2)./m);
% x4=(s4+x2(1,1))/2;



   X=[x1,x2,y1,x3];
   Y=[y1,x2,x1,y3];
   
   
   scatter(X,Y,8,'ko');
   axis equal;