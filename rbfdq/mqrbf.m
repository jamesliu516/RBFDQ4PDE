function [rder]=mqrbf(pxy,xy,c)
% This program is used to calculate the derivative coefficients in the Local MQ-DQ method.
% c---- INPUT: pxy, xy, c , xy is a  row vector like (0.3,0.5)
% c---- OUTPUT: r
% c---- pxy: store the positions of the supporting points np*2 matrix
% c---- xy: store the position of the reference node
% c---- c: shape parameter for the MQ radial basis function 
% c---- rder: vector of computed derivative coefficients
% c---- Some important symbols and variables
% c---- np: the number of supporting points
% c---- A: coefficient matrix constructed from the basis functions
% c---- b: derivative vectors of the basis functions 

np=size(pxy,1);
nd=np+1;
rder=zeros(nd,5);
pn=zeros(nd,2);

for i=1:nd
    if i~= nd
        pn(i,:)=pxy(i,:);
    else
        pn(i,:)=xy;
    end
end

scaling=0.0;

for i=1:np
    dx=pxy(i,1)-xy(1);
    dy=pxy(i,2)-xy(2);
    scaling = max(scaling,sqrt(dx^2+dy^2));
end

scaling = scaling*2;

a=zeros(nd,nd);
b=zeros(nd,5);

a(nd,:)=1.0;

for i=1:nd-1
    for j=1:nd
        dx=(pn(j,1)-pn(i,1))/scaling;       
        dy=(pn(j,2)-pn(i,2))/scaling;
        dxk=(pn(j,1)-pn(nd,1))/scaling;
        dyk=(pn(j,2)-pn(nd,2))/scaling;
        a(i,j)=sqrt(dx*dx+dy*dy+c)-sqrt(dxk*dxk+dyk*dyk+c);
    end
end

for i=1:nd-1
    dx=(-pn(i,1)+pn(nd,1))/scaling;
    dy=(-pn(i,2)+pn(nd,2))/scaling;
    ffunc=sqrt(dx*dx+dy*dy+c);
    b(i,1)=dx/ffunc;
    
    b(i,2)=dy/ffunc;
    b(i,3)=(dy*dy+c)/(ffunc^3.0)-1.d0/sqrt(c);
    b(i,5)=-dx*dy/(ffunc^3.0);
    b(i,4)=(dx*dx+c)/(ffunc^3.)-1.d0/sqrt(c);
end

b(nd,1)=0.0;
b(nd,2)=0.0;
b(nd,3)=0.0;
b(nd,4)=0.0;
b(nd,5)=0.0;

for i=1:nd
    for j=1:5
        rder(i,j)=b(i,j);
    end
end

rder=a\rder;


for ik1=1:5
    for ik2=1:nd
        if (ik1==1 || ik1==2)
            rder(ik2,ik1)=rder(ik2,ik1)/scaling;            
        elseif (ik1==3 || ik1==4 || ik1==5)
            rder(ik2,ik1)=rder(ik2,ik1)/scaling/scaling;
        end
    end
end








        
