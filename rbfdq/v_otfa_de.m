%variable-order time fractional advection-diffusion equation solver
%v_otfa_de
clear 
clc
close all
hold off
% test mqrbf and meshfree grid treatment
%
global ppp meshden  pointboun typPoints domain racLow racHigh
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
meshden=0.1;
thet=1;  % theta method


npoin=size(ppp,1);

examp=1; %different case
domain=1; %1 [-1,1]*[-1,1],2: unit circle
if domain==1
    racLow=[0,0]; % left down
    racHigh=[1,1];  % right up
end
    
boundType=1; %1 Dirichlet, 2 Neumann

NtimeStep=100;
Tend=1.0;
dlt=Tend/NtimeStep;
Tnow=0;

unum=zeros(npoin,NtimeStep+1);
Fnum=zeros(npoin,1);

switch examp
    case 1
        kapaFun=@(x,y,t) (1  );
        vecSp1=@(x,y,t) (1 );
        vecSp2=@(x,y,t) (1 );        
        vo_alpha=@(x,y,t) (0.8-0.1*cos(x*t).* sin(x)-0.1*cos(y*t).*sin(y));
      %  vo_alpha=@(x,y,t) (0.5);
        sourceF=@(x,y,t) (2*t.^(2-vo_alpha(x,y,t))./gamma(3-vo_alpha(x,y,t)) ...
            +2*x+2*y-4);
        uexact=@(x,y,t) (x.^2+y.^2+t.^2);
    otherwise
        warning('Unexpected example.');
        pause;
end

muFun=@(x,y,t) (dlt.^ vo_alpha(x,y,t).* gamma(2-vo_alpha(x,y,t)));
bb=@(x,y,t,jj) ((jj+1).^(1-vo_alpha(x,y,t))-jj.^(1-vo_alpha(x,y,t)));

meshfreeTreat;

pxy=cell(npoin,1);
for ipoin=1:npoin
    for jk=1:n_pointPoint2(ipoin)
       pxy{ipoin}=[pxy{ipoin}; ppp(pointsPoint2(ipoin,jk),:)];
    end
end


rder=cell(npoin,1);
c=5;
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
% c---- rder(:,1) coefficients for u_x, rder(:,2) for u_y, 
% c---- rder(:,3) d^2u/dx^2  rder(:,4) d^2u/dy^2   rder(:,5) d^2u/dxdy


nStep=0;
%initial data from time =0
% for ipoin=1:npoin
%     unum(ipoin,1)=uexact(ppp(ipoin,1), ppp(ipoin,2), 0);
% end
unum(:,1)=uexact(ppp(:,1), ppp(:,2), 0);

while Tnow<Tend 
    acoe=zeros(npoin,npoin);
    Tnow=Tnow+dlt;
    nStep=nStep+1;
%     if Tnow>Tend
%         break
%     end
    
    for ipoin=1:npoin
        if typPoints(ipoin)==0
            att1=rder{ipoin};
            acoe(ipoin,ipoin)=1;
            for jk=1:n_pointPoint2(ipoin)
                nbpoin=pointsPoint2(ipoin,jk);
                rt=-(att1(jk,3)+att1(jk,4)) ...
                    *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
                acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                    +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow);                              
            end
            
            % att1(n_pointPoint2(ipoin)+1,1)
            jk=n_pointPoint2(ipoin)+1;
            
            nbpoin=ipoin;
            rt=-(att1(jk,3)+att1(jk,4)) ...
                *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
            acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow);
            tmpsum=0;
            for itp=1:nStep-2
                tmpsum=tmpsum+(bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,itp) ...
                    -bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,itp+1))*unum(ipoin,nStep-1-itp+1);
            end
            
            Fnum(ipoin)=(1-bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,1)) ...
                *unum(ipoin,nStep-1+1)+ tmpsum+bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep) ...
                *unum(ipoin,0+1)+muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                *sourceF(ppp(ipoin,1),ppp(ipoin,2),Tnow);
        end
        
        if  typPoints(ipoin)==1
            acoe(ipoin,ipoin)=1;
            Fnum(ipoin)=uexact(ppp(ipoin,1),ppp(ipoin,2),Tnow);
        end                                               
    end   
    unum(:,nStep+1)=acoe\Fnum;
end
        
        
        
uerr=unum(:,NtimeStep+1)-uexact(ppp(:,1),ppp(:,2),Tend);  
rmserr=sqrt(sum(uerr.^2)/npoin);
        
    


figure(2)
hold off

plot(ppp(:,1),ppp(:,2),'b.','MarkerSize',20);
hold on
plot(ppp(pointboun,1),ppp(pointboun,2),'r.','MarkerSize',20);
xlabel('x'); ylabel('y');

if domain==1
    axis([racLow(1),racHigh(1),racLow(2),racHigh(2)])
end
axis equal

figure(3)
plot3(ppp(:,1),ppp(:,2), unum(:,NtimeStep+1), 'b.','MarkerSize',20);
xlabel('x'); ylabel('y');
zlabel('u^h({\bf x}, T)');
grid on

figure(4)
plot3(ppp(:,1),ppp(:,2), abs(uerr), 'b.','MarkerSize',20);
zlabel('Error');
xlabel('x'); ylabel('y');
grid on
