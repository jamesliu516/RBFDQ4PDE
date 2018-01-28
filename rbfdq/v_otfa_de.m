%variable-order time fractional advection-diffusion equation solver
%v_otfa_de
clear 
clc
%close all
hold off
% test mqrbf and meshfree grid treatment
%
global ppp meshden  pointboun typPoints domain racLow racHigh
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
meshden=0.1;

examp=1; %different case
domain=2; %1 [0,1]*[0,1],2: unit circle
if domain==1
    racLow=[0,0]; % left down
    racHigh=[1,1];  % right up
end
boundType=2; %1 Dirichlet, 2 Neumann


cellBool=0; % 1: cell 0: map
HRBFDQ=0; %0: rbf dq by Shu, 1: hermite RBFDQ
c=25; 
meshfreeTreat;

thet=1;  % theta method

npoin=size(ppp,1);
                                 
NtimeStep=200;
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
        vo_alpha=@(x,y,t) (0.8-0.1*cos(x.*t).* sin(x)-0.1*cos(y.*t).*sin(y));
        %vo_alpha=@(x,y,t) (0.8);
        sourceF=@(x,y,t) (2*t.^(2-vo_alpha(x,y,t))./gamma(3-vo_alpha(x,y,t)) ...
            +2*x+2*y-4);
        uexact=@(x,y,t) (x.^2+y.^2+t.^2);
        dfx1=@(x,y,t) (2*x);
        dfy1=@(x,y,t) (2*y);              
    otherwise
        warning('Unexpected example.');
        pause;
end

muFun=@(x,y,t) (dlt.^ vo_alpha(x,y,t).* gamma(2-vo_alpha(x,y,t)));
bb=@(x,y,t,jj) ((jj+1).^(1-vo_alpha(x,y,t))-jj.^(1-vo_alpha(x,y,t)));


typPoints(pointboun)=boundType; %%all boundary points: Neumann boundary points

pxy=cell(npoin,1);
for ipoin=1:npoin
    for jk=1:n_pointPoint2(ipoin)
       pxy{ipoin}=[pxy{ipoin}; ppp(pointsPoint2(ipoin,jk),:)];
    end
end

numbp=size(pointboun,1);


if cellBool==0
    rdernb=cell(numbp,1);
    rdernb1=cell(numbp,1);
    rdernb2=cell(numbp,1);
    rdernbMap=containers.Map(pointboun,rdernb);
    pxynbMap=containers.Map(pointboun,rdernb1);
    pxynbnorMap=containers.Map(pointboun,rdernb2);
end

if cellBool==1
    rdernbMap=cell(npoin,1);
    pxynbMap=cell(npoin,1);
    pxynbnorMap=cell(npoin,1);
end
           

rder=cell(npoin,1);


for ipoin=1:npoin    
    pxy11=pxy{ipoin};
    xy=ppp(ipoin,:);
    rd=mqrbf(pxy11,xy,c);
    
    rder{ipoin}=[rder{ipoin}; rd];
    if  typPoints(ipoin)==2
        pnxy1=pointsPoint2(ipoin,1:n_pointPoint2(ipoin));
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
        if cellBool==0
            if isKey(rdernbMap,ipoin)
                rdernbMap(ipoin)=[rdernbMap(ipoin);rd2];
                pxynbMap(ipoin)=[pxynbMap(ipoin);pxynb];
                pxynbnorMap(ipoin)=[pxynbnorMap(ipoin);pxynbnor];
            else
                disp('Wrone rdernbMap');
                pause;
            end           
        end
        
        if cellBool==1
            rdernbMap{ipoin}=[rdernbMap{ipoin};rd2];
            pxynbMap{ipoin}=[pxynbMap{ipoin};pxynb];
            pxynbnorMap{ipoin}=[pxynbnorMap{ipoin};pxynbnor];
        end
            
                 
    end
    
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
%%initial data from time =0
% for ipoin=1:npoin
%     unum(ipoin,1)=uexact(ppp(ipoin,1), ppp(ipoin,2), 0);
% end
unum(:,1)=uexact(ppp(:,1), ppp(:,2), 0);

while Tnow<Tend 
    acoe=zeros(npoin,npoin);
    Fnum=zeros(npoin,1);
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
            %nStep == k+1
            for itp=1:nStep-2
                tmpsum=tmpsum+(bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,itp) ...
                    -bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,itp+1))*unum(ipoin,nStep-1-itp+1);
            end
            % nStep-1 ?????
            Fnum(ipoin)=(1-bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,1)) ...
                *unum(ipoin,nStep-1+1)+ tmpsum+bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep) ...
                *unum(ipoin,0+1)+muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                *sourceF(ppp(ipoin,1),ppp(ipoin,2),Tnow);
        end
        
        if  typPoints(ipoin)==1          
            acoe(ipoin,ipoin)=1;
            Fnum(ipoin)=uexact(ppp(ipoin,1),ppp(ipoin,2),Tnow);           
        end
        %%%neumann boundary can only do for unit circle  
        if  typPoints(ipoin)==2
            if cellBool==0
                rd2=rdernbMap(ipoin);
                pxynb=pxynbMap(ipoin);
                pxynbnor=pxynbnorMap(ipoin);
            end
            
            if cellBool==1
                rd2=rdernbMap{ipoin};
                pxynb=pxynbMap{ipoin};
                pxynbnor=pxynbnorMap{ipoin};
            end
            
            nor=ppp(ipoin,:);
            
            if HRBFDQ==1
                % rt=0.0;
                for jk=1:n_pointPoint2(ipoin)
                    %rt=rt+(rd2(jk,1)*nor(1)+rd2(jk,2)*nor(2))*af(pointsPoint2(ipoin,jk));
                    acoe(ipoin,pointsPoint2(ipoin,jk))=rd2(jk,1)*nor(1)+rd2(jk,2)*nor(2);
                end
                % xy1=ppp(ipoin,:);
                nd=n_pointPoint2(ipoin)+1;
                acoe(ipoin,ipoin)=rd2(nd,1)*nor(1)+rd2(nd,2)*nor(2);
                % npnb=size(rd2,1)-nd;
                
                npnb=size(pxynb,1);
                rt=0;
                for jk=1:npnb
                    tm1=dfx1(pxynb(jk,1),pxynb(jk,2),Tnow)*pxynbnor(jk,1);
                    tm2=dfy1(pxynb(jk,1),pxynb(jk,2),Tnow)*pxynbnor(jk,2);
                    rt=rt+(rd2(jk+nd,1)*nor(1)+rd2(jk+nd,2)*nor(2))*(tm1+tm2);
                end
                xy1=ppp(ipoin,:);
                tm1= dfx1(xy1(1),xy1(2),Tnow)*nor(1)+dfy1(xy1(1),xy1(2),Tnow)*nor(2);
                
                Fnum(ipoin)=tm1-rt;
            end
            
            if HRBFDQ==0
                rd3=rder{ipoin};
                for jk=1:n_pointPoint2(ipoin)
                    acoe(ipoin,pointsPoint2(ipoin,jk))= ...
                        rd3(jk,1)*nor(1)+rd3(jk,2)*nor(2);
                end
                nd=n_pointPoint2(ipoin)+1;
                
                acoe(ipoin,ipoin)=rd3(nd,1)*nor(1)+rd3(nd,2)*nor(2);
                xy1=ppp(ipoin,:);
                tm1= dfx1(xy1(1),xy1(2),Tnow)*nor(1)+dfy1(xy1(1),xy1(2),Tnow)*nor(2);
                
                Fnum(ipoin)=tm1;
                
            end                                      
        end
    end
    unum(:,nStep+1)=acoe\Fnum;
end
        
        
        
uerr=unum(:,NtimeStep+1)-uexact(ppp(:,1),ppp(:,2),Tend);  
rmserr=sqrt(sum(uerr.^2)/npoin);
lwqerr=max(abs(uerr));
l2err=sqrt(sum(uerr.^2)/sum(uexact(ppp(:,1),ppp(:,2),Tend).^2));
        
    


figure(2)
hold off

plot(ppp(:,1),ppp(:,2),'b.','MarkerSize',20);
hold on
plot(ppp(pointboun,1),ppp(pointboun,2),'r.','MarkerSize',20);
xlabel('x'); ylabel('y');

if domain==1
    axis([racLow(1),racHigh(1),racLow(2),racHigh(2)])
end
if domain ==2
    axis equal
end

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
