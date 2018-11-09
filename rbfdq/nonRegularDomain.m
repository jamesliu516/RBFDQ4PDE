%nonregularDomain
%plot a non regular 

the1=0:0.0125:2*pi;
npts=6;
r=1.0/npts/npts *(1+2*npts+npts^2-(npts+1)*cos(npts*the1));

x=r.*cos(the1);

y=r.*sin(the1);

plot(x,y)
ze=zeros(length(x),1);
xy=[x' y' ze];
