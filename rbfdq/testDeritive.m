clear
clc
% test mqrbf and meshfree grid treatment
global ppp meshden  pointboun typPoints
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
meshden=0.05;

meshfreeTreat;

% f=@(x,y) (x.^2-x.*y.^3+4);
% 
% dfx1=@(x,y) (2*x-y.^3);
% dfy1=@(x,y) (-3*x.*y.^2);
% 
% f=@(x,y) (y.^2+y.*x.^3+4);
% 
% dfy1=@(x,y) (2*y+x.^3);
% dfx1=@(x,y) (3*y.*x.^2);


f=@(x,y) (x.^3+y.^3+x+y+6);
dfx1=@(x,y) (3*x.^2+1);
dfy1=@(x,y) (3*y.^2+1);

npoin=size(ppp,1);
af=f(ppp(:,1),ppp(:,2));

adfx1=dfx1(ppp(:,1),ppp(:,2));
adfy1=dfy1(ppp(:,1),ppp(:,2));
cadfy1=zeros(npoin,1);
cadfx1=zeros(npoin,1);

pxy=cell(npoin,1);
for ipoin=1:npoin
    for jk=1:n_pointPoint2(ipoin)
       pxy{ipoin}=[pxy{ipoin}; ppp(pointsPoint2(ipoin,jk),:)];
    end
end


rder=cell(npoin,1);
c=25;
for ipoin=1:npoin    
    pxy11=pxy{ipoin};
    xy=ppp(ipoin,:);
    rd=mqrbf(pxy11,xy,c);
    rder{ipoin}=[rder{ipoin}; rd];
end
%rt=0;
sumerr=0.0;
for ipoin=1:npoin
    att=rder{ipoin};
    rt=0.0;
    for jk=1:n_pointPoint2(ipoin)
       rt=rt+ att(jk,1)*af(pointsPoint2(ipoin,jk));
    end
    rt=rt+att(n_pointPoint2(ipoin)+1,1)*af(ipoin);
    sumerr=sumerr+((rt-adfx1(ipoin))/(abs(adfx1(ipoin))+1e-8))^2;
end

sumerr=sqrt(sumerr/npoin);

rt=0;
sumerr2=0.0;
for ipoin=1:npoin
    att1=rder{ipoin};
    rt=0.0;
    for jk=1:n_pointPoint2(ipoin)
       rt=rt+ att1(jk,2)*af(pointsPoint2(ipoin,jk));
    end
    rt=rt+att1(n_pointPoint2(ipoin)+1,2)*af(ipoin);
    cadfy1(ipoin)=rt;
    sumerr2=sumerr2+((rt-adfy1(ipoin))/(abs(adfy1(ipoin))+1e-8))^2;  
end

 sumerr2=sqrt(sumerr2/npoin);
 
% errvec=abs(cadfy1-adfy1);
% [err,id]=max(errvec); 
% 
% sum(pointboun==id)
% find(pointboun==id)
% 
% max(errvec)
%     
    
        
        



%function [rder]=mqrbf(pxy,xy,c)





