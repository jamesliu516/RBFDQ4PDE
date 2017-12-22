function [rder,rderN]=mqrbfNB(pxy,pxynb, pxynbnor, xy, c)
% This program is used to calculate the derivative coefficients in the Local MQ-DQ method.
% c---- INPUT: pxy, xy, c , xy is a  row vector like (0.3,0.5)
% c----- pxy do not include point xy
% c---- OUTPUT: r
% c---- pxy: store the positions of the supporting points np*2 matrix
% c---- xy: store the position of the reference node
% c---- pxynb Neumann boundary points at the boundary
% c---- pxynbnor Normal direction at Neumann boundary points
% c---- c: shape parameter for the MQ radial basis function 
% c---- rder: vector of computed derivative coefficients
% c---- Some important symbols and variables
% c---- np: the number of supporting points
% c---- A: coefficient matrix constructed from the basis functions
% c---- b: derivative vectors of the basis functions 

np=size(pxy,1);
npnb=size(pxynb,1);
nd=np+1;
nd1=np+1+npnb;
rder=zeros(nd,5);
pn=zeros(nd,2);

for ji=1:nd
    if ji~= nd
        pn(ji,:)=pxy(ji,:);
    else
        pn(ji,:)=xy;
    end
end


scaling=0.0;

for ji=1:np
    dx=pxy(ji,1)-xy(1);
    dy=pxy(ji,2)-xy(2);
    scaling = max(scaling,sqrt(dx^2+dy^2));
end

scaling = scaling*2.0;

a=zeros(nd1,nd1);
b=zeros(nd1,2);

a(nd,:)=1.0;

for ii=1:nd-1
    for jj=1:nd
        dx=(pn(jj,1)-pn(ii,1))/scaling;       
        dy=(pn(jj,2)-pn(ii,2))/scaling;
        dxk=(pn(jj,1)-pn(nd,1))/scaling;
        dyk=(pn(jj,2)-pn(nd,2))/scaling;
        a(ii,jj)=sqrt(dx*dx+dy*dy+c)-sqrt(dxk*dxk+dyk*dyk+c);
    end
end

for ii=1:nd-1
    dx=(-pn(ii,1)+pn(nd,1))/scaling;
    dy=(-pn(ii,2)+pn(nd,2))/scaling;
    ffunc=sqrt(dx*dx+dy*dy+c);
    b(ii,1)=dx/ffunc;
    
    b(ii,2)=dy/ffunc;
    b(ii,3)=(dy*dy+c)/(ffunc^3.0)-1.d0/sqrt(c);
    b(ii,5)=-dx*dy/(ffunc^3.0);
    b(ii,4)=(dx*dx+c)/(ffunc^3.)-1.d0/sqrt(c);
end

b(nd,1)=0.0;
b(nd,2)=0.0;
b(nd,3)=0.0;
b(nd,4)=0.0;
b(nd,5)=0.0;

for ii=1:nd
    for jj=1:5
        rder(ii,jj)=b(ii,jj);
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








        
