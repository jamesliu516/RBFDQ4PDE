function [rder]=mqrbf(pxy,xy,c)
% This program is used to calculate the derivative coefficients in the Local MQ-DQ method.
% c---- INPUT: pxy, xy, c , xy is a  row vector like (0.3,0.5)
% c---- OUTPUT: r
% c---- pxy: store the positions of the supporting points np*2 matrix
% c---- xy: store the position of the reference node
% c---- c: shape parameter for the MQ radial basis function 
% c---- r: vector of computed derivative coefficients
% c---- Some important symbols and variables
% c---- np: the number of supporting points
% c---- A: coefficient matrix constructed from the basis functions
% c---- b: derivative vectors of the basis functions 

np=size(pxy,1);
nd=np+1;
rder=zeros(nd,5);
pn=zeros(nd,2);

a=zeros(nd,nd);
b=zeros(nd,5);

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

scaling = scaling *2;

a=1;



        
