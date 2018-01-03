clear
clc
clf 
hold off
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
sumexact=0.0;
for ipoin=1:npoin
    att=rder{ipoin};
    rt=0.0;
    for jk=1:n_pointPoint2(ipoin)
       rt=rt+ att(jk,1)*af(pointsPoint2(ipoin,jk));
    end
    rt=rt+att(n_pointPoint2(ipoin)+1,1)*af(ipoin);
    cadfx1(ipoin)=rt;    
    sumerr=sumerr+(rt-adfx1(ipoin))^2;
    sumexact=sumexact+abs(adfx1(ipoin))^2;
end

sumerr=sqrt(sumerr/sumexact);

rt=0;
sumerr2=0.0;
sumexact2=0.0;
for ipoin=1:npoin
    att1=rder{ipoin};
    rt=0.0;
    for jk=1:n_pointPoint2(ipoin)
       rt=rt+ att1(jk,2)*af(pointsPoint2(ipoin,jk));
    end
    rt=rt+att1(n_pointPoint2(ipoin)+1,2)*af(ipoin);
    cadfy1(ipoin)=rt;
    sumerr2=sumerr2+(rt-adfy1(ipoin))^2; 
    sumexact2=sumexact2+abs(adfy1(ipoin))^2;  
end

 sumerr2=sqrt(sumerr2/sumexact2);
 
% errvec=abs(cadfy1-adfy1);
% [err,id]=max(errvec); 
% 
% sum(pointboun==id)
% find(pointboun==id)
% 
% max(errvec)
%     
    
figure(1)

% for ipoin=1:npoin
%     if typPoints(ipoin)==1
%         pcol='r.';
%     else
%         pcol='b.';
%     end
    plot(ppp(:,1),ppp(:,2),'b.','MarkerSize',20);
    hold on
    plot(ppp(pointboun,1),ppp(pointboun,2),'r.','MarkerSize',20);
    
% end
axis equal

% test Neumann boundary condition
%function [rder]=mqrbfNB(pxy,pxynb, pxynbnor, xy, c)

sumerr3=0.0;
sumdudn=0.0;
pnberr=[];
for k=1:size(pointboun,1)
    ipoin=pointboun(k);
    
    pnxy1=pointsPoint2(ipoin,1:n_pointPoint2(ipoin));
%     if length(pnxy1)~=n_pointPoint2(ipoin)
%         disp('error pnxy1');
%     end
    pxy1=ppp(pnxy1',:);
    xy1=ppp(ipoin,:);
    pxynb=[];
    pxynbnor=[];
    pxynb=[pxynb; xy1];
    %nor=(xy1-0.0)/norm(xy1);
    nor=xy1;
    pxynbnor=[pxynbnor; nor];
    for m=1:n_pointPoint2(ipoin)
       pnb =find(pointboun==pointsPoint2(ipoin,m));
       if ~isempty(pnb)
           pxynb=[pxynb; ppp(pointboun(pnb),:)];
           xytm=ppp(pointboun(pnb),:); % here circle r =1
           pxynbnor=[pxynbnor;xytm];
       end
    end
    
    rd2= mqrbfNB(pxy1,pxynb, pxynbnor, xy1, c);
    rt=0.0;
    
    for jk=1:n_pointPoint2(ipoin)
        rt=rt+(rd2(jk,1)*nor(1)+rd2(jk,2)*nor(2))*af(pointsPoint2(ipoin,jk));
    end
    nd=n_pointPoint2(ipoin)+1;
    
    rt=rt+(rd2(nd,1)*nor(1)+rd2(nd,2)*nor(2))*af(ipoin);
    npnb=size(pxynb,1);

    for jk=1:npnb
        tm1=dfx1(pxynb(jk,1),pxynb(jk,2))*pxynbnor(jk,1);
        tm2=dfy1(pxynb(jk,1),pxynb(jk,2))*pxynbnor(jk,2); 
        rt=rt+(rd2(jk+nd,1)*nor(1)+rd2(jk+nd,2)*nor(2))*(tm1+tm2);
    end
    tm1= dfx1(xy1(1),xy1(2))*nor(1)+dfy1(xy1(1),xy1(2))*nor(2);
  %  rt-tm1
    
    
    sumerr3=sumerr3+(rt-tm1)^2;
    sumdudn=sumdudn+tm1^2;
    pnberr=[pnberr; rt-tm1];
end

sumerr3=sqrt(sumerr3/sumdudn);


        
        
    
    
    
    
        
        
    
   
%    
%   sumerr=0.0;
% sumexact=0.0;
% for ipoin=1:npoin
%     att=rder{ipoin};
%     rt=0.0;
%     for jk=1:n_pointPoint2(ipoin)
%        rt=rt+ att(jk,1)*af(pointsPoint2(ipoin,jk));
%     end
%     rt=rt+att(n_pointPoint2(ipoin)+1,1)*af(ipoin);
%     cadfx1(ipoin)=rt;    
%     sumerr=sumerr+(rt-adfx1(ipoin))^2;
%     sumexact=sumexact+abs(adfx1(ipoin))^2;
% end
% 
% sumerr=sqrt(sumerr/sumexact);  
    
        
    
    
    
    




