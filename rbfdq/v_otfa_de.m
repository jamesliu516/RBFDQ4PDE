%variable-order time fractional advection-diffusion equation solver
%v_otfa_de
clear
clc
clf 
hold off
% test mqrbf and meshfree grid treatment
global ppp meshden  pointboun typPoints domain racLow racHigh
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
meshden=0.05;
thet=1;  % theta method


npoin=size(ppp,1);

examp=1; %different case
domain=1; %1 [-1,1]*[-1,1],2: unit circle
if domain==1
    racLow=[0,0]; % left down
    racHigh=[1,1];  % right up
end
    
boundType=1; %1 Dirichlet, 2 Neumann
dlt=0.01;

switch examp
    case 1
        kapaFun=@(x,y,t) (1  );
        vecSp=@(x,y,t) (1 );
        vo_alpha=@(x,y,t) (0.8-0.1*cos(x*t).* sin(x)-0.1*cos(y*t).*sin(y));
        sourceF=@(x,y,t) (2*t.^(2-vo_alpha(x,y,t))./gamma(3-vo_alpha(x,y,t)) ...
            +2*x+2*y-4);
        uexact=@(x,y,t) (x.^2+y.^2+t.^2);
    otherwise
        warning('Unexpected example.');
        pause;
end

muFun=@(x,y,t) (dlt.^ vo_alpha(x,y,t).* gamma(2-vo_alpha(x,y,t)));

meshfreeTreat;

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

% for ipoin=1:npoin
%     att1=rder{ipoin};
%     rt=0.0;
%     for jk=1:n_pointPoint2(ipoin)
%        rt=rt+ att1(jk,2)*af(pointsPoint2(ipoin,jk));
%     end
%     rt=rt+att1(n_pointPoint2(ipoin)+1,2)*af(ipoin);
%     cadfy1(ipoin)=rt;
%     sumerr2=sumerr2+(rt-adfy1(ipoin))^2; 
%     sumexact2=sumexact2+abs(adfy1(ipoin))^2;  
% end

acoe=zeros(npoin,npoin);
for ipoin=1:npoin
    if typPoints(ipoin)==0
        att1=rder{ipoin};
        for jk=1:n_pointPoint2(ipoin)
          %  rt=rt+ att1(jk,2)*af(pointsPoint2(ipoin,jk));
          
          
        end
        
        
        
    


figure(2)
hold off

plot(ppp(:,1),ppp(:,2),'b.','MarkerSize',20);
hold on
plot(ppp(pointboun,1),ppp(pointboun,2),'r.','MarkerSize',20);

if domain==1
    axis([racLow(1),racHigh(1),racLow(2),racHigh(2)])
end

