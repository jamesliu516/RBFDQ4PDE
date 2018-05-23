%variable-order time fractional advection-diffusion equation solver
%v_otfa_de
clear 
clc
%close all
hold off
% test mqrbf and meshfree grid treatment
%
global ppp meshden  pointboun typPoints domain racLow racHigh ttt
%pointboun: boundary node number

%global n_pointPoint pointsPoint

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
HRBFDQ=1; %0: rbf dq by Shu, 1: hermite RBFDQ
%boundInEq=0; % 1 include boundary point Eq, 0 no

c=5; 

global su2mesh
su2mesh=1;
meshfreeTreat;

thet=1;  % theta method

npoin=size(ppp,1);
zeroORnpoin=0; %% if boundary points are included in eqs least square method is used
                                             
NtimeStep=200;
Tend=1;
dlt=Tend/NtimeStep;
Tnow=0;



switch examp
    case 1
        kapaFun=@(x,y,t) (1  );
        vecSp1=@(x,y,t) (1 );
        vecSp2=@(x,y,t) (1 );        
        vo_alpha=@(x,y,t) (0.8-0.1*cos(x.*t).* sin(x)-0.1*cos(y.*t).*sin(y));
       % vo_alpha=@(x,y,t) (1.0);
        sourceF=@(x,y,t) (2*t.^(2-vo_alpha(x,y,t))./gamma(3-vo_alpha(x,y,t)) ...
            +2*x+2*y-4);
        uexact=@(x,y,t) (x.^2+y.^2+t.^2);
        dfx1=@(x,y,t) (2*x);
        dfy1=@(x,y,t) (2*y);              
    otherwise
        warning('Unexpected example.');
        pause;
end

muFun=@(x,y,t,dlt) (dlt.^ vo_alpha(x,y,t).* gamma(2-vo_alpha(x,y,t)));
bb=@(x,y,t,jj) ((jj+1).^(1-vo_alpha(x,y,t))-jj.^(1-vo_alpha(x,y,t)));


typPoints(pointboun)=boundType; %%all boundary points: Neumann boundary points

pxy=cell(npoin,1);
for ipoin=1:npoin
    for jk=1:n_pointPoint2(ipoin)
       pxy{ipoin}=[pxy{ipoin}; ppp(pointsPoint2(ipoin,jk),:)];
    end
end
numbp=size(pointboun,1);
nmlPboun=zeros(numbp,2); %normal direction over the boundary points
g2locNmlPbnd=containers.Map(pointboun,zeros(numbp,1));
neumBndV=zeros(numbp,1);
unum=zeros(npoin+numbp,NtimeStep+1);
Fnum=zeros(npoin+numbp,1);

for ipb=1:numbp
    nmlPboun(ipb,:)=ppp(pointboun(ipb),:)+1e-6; % only right for unit circle and origin is the center 
    g2locNmlPbnd(pointboun(ipb))=ipb;
end


tmpCell=cell(npoin,3);   %tmpCell(0(1), [... number on boundary ],  [... mqrbfNb coefficents]
for ipoin=1:npoin
    tmpCell{ipoin,1}=0;  
    tmpCell{ipoin,2}=[];
    tmpCell{ipoin,3}=[];
    if typPoints(ipoin)==2
        tmpCell{ipoin,1}=1; 
        tmpCell{ipoin,2}=[tmpCell{ipoin,2}, ipoin];
    end
    
    for m=1:n_pointPoint2(ipoin)
        pnb =find(pointboun==pointsPoint2(ipoin,m));
        if ~isempty(pnb)
            tmpCell{ipoin,1}=1;
            tmpCell{ipoin,2}=[tmpCell{ipoin,2},pointboun(pnb)];   
        end
    end    
end

% rd2= mqrbfNB(pxy1,xy1,pxynb, pxynbnor, c);
for ipoin=1:npoin
    pxy11=pxy{ipoin};
    xy=ppp(ipoin,:);
    pxynb=[];
    pxynbnor=[];
    if tmpCell{ipoin,1}==1 % 1 denotes the present point or some neighbor point on the Neumann boundary
        nbpArray=tmpCell{ipoin,2};
        for itm=1:size(nbpArray,2)
            pxynb=[pxynb; ppp(nbpArray(itm),:)];
            pxynbnor=[pxynbnor; nmlPboun(g2locNmlPbnd(nbpArray(itm)),:)];
        end
        
        rdrd= mqrbfNB(pxy11,xy,pxynb, pxynbnor, c);
    elseif tmpCell{ipoin,1}==0
        rdrd= mqrbf(pxy11,xy, c);          
    end
    tmpCell{ipoin,3}=[tmpCell{ipoin,3}, rdrd];   
end

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
        nnbt=1;
        % only one Neumann boundary point then comment following for cycle
        for m=1:n_pointPoint2(ipoin)
            pnb =find(pointboun==pointsPoint2(ipoin,m));
            if ~isempty(pnb)
                pxynb=[pxynb; ppp(pointboun(pnb),:)];
                xytm=ppp(pointboun(pnb),:); % here circle r =1
                pxynbnor=[pxynbnor;xytm];
                nnbt=nnbt+1;
            end
%             if nnbt >3
%                 break;
%             end
        end
        
        rd2= mqrbfNB(pxy1,xy1,pxynb, pxynbnor, c);
        if cellBool==0
            if isKey(rdernbMap,ipoin)
                rdernbMap(ipoin)=[rdernbMap(ipoin);rd2];
                pxynbMap(ipoin)=[pxynbMap(ipoin);pxynb];
                pxynbnorMap(ipoin)=[pxynbnorMap(ipoin);pxynbnor];
            else
                disp('Wrong rdernbMap');
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
unum(1:npoin,1)=uexact(ppp(:,1), ppp(:,2), 0);

while Tnow<Tend
    acoe=zeros(npoin+numbp,npoin+numbp);
    Fnum=zeros(npoin+numbp,1);
    Tnow=Tnow+dlt;
    nStep=nStep+1;
    %     if Tnow>Tend
    %         break
    %     end
    
   % nBoundEq=0;
    for ibn=1:numbp
        xy1=ppp(pointboun(ibn),:);
        nor=nmlPboun(ibn,:);
        neumBndV(ibn)=dfx1(xy1(1),xy1(2),Tnow)*nor(1) ...
            +dfy1(xy1(1),xy1(2),Tnow)*nor(2);
    end
    
    for ipoin=1:npoin
        %if typPoints(ipoin)==0 || (typPoints(ipoin)==2 && boundInEq==1)
        if tmpCell{ipoin,1}==0
            att1=rder{ipoin};
            acoe(ipoin,ipoin)=1;
            for jk=1:n_pointPoint2(ipoin)
                nbpoin=pointsPoint2(ipoin,jk);
                rt=-(att1(jk,3)+att1(jk,4)) ...
                    *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
                acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                    +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt);
            end
            
            % att1(n_pointPoint2(ipoin)+1,1)
            jk=n_pointPoint2(ipoin)+1;
            
            nbpoin=ipoin;
            rt=-(att1(jk,3)+att1(jk,4)) ...
                *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
            acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt);
            tmpsum=0;

            for itp=1:nStep-1
                tmpsum=tmpsum+(bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-itp-1) ...
                    -bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-itp)) ...
                    *unum(ipoin,itp+1);
            end   
            Fnum(ipoin)=bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-1) ...
                * unum(ipoin,1)+tmpsum+muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt) ...
                 *sourceF(ppp(ipoin,1),ppp(ipoin,2),Tnow);

        end
               
        if  typPoints(ipoin)==1
            acoe(ipoin,ipoin)=1;
            Fnum(ipoin)=uexact(ppp(ipoin,1),ppp(ipoin,2),Tnow);
        end
        %%%neumann boundary can only do for unit circle  
      %  if tmpCell{ipoin,1}==1 && typPoints(ipoin)~=2
        if tmpCell{ipoin,1}==1           
            att1=tmpCell{ipoin,3};
            acoe(ipoin,ipoin)=1;
            for jk=1:n_pointPoint2(ipoin)
                nbpoin=pointsPoint2(ipoin,jk);
                rt=-(att1(jk,3)+att1(jk,4)) ...
                    *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
                acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                    +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt);
            end
            
            % att1(n_pointPoint2(ipoin)+1,1)
            jk=n_pointPoint2(ipoin)+1;
            
            nbpoin=ipoin;
            rt=-(att1(jk,3)+att1(jk,4)) ...
                *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
            acoe(ipoin,nbpoin)=acoe(ipoin,nbpoin) ...
                +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt);
            tmpsum=0; 
            
            for itp=1:nStep-1
                tmpsum=tmpsum+(bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-itp-1) ...
                    -bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-itp)) ...
                    *unum(ipoin,itp+1);
            end   
            Fnum(ipoin)=bb(ppp(ipoin,1),ppp(ipoin,2),Tnow,nStep-1) ...
                * unum(ipoin,1)+tmpsum+muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt) ...
                 *sourceF(ppp(ipoin,1),ppp(ipoin,2),Tnow); 
             
            npsp=n_pointPoint2(ipoin)+1; 
            bndp=tmpCell{ipoin,2};
         %   bndnmlVadd=0.0;
            %neumBndV(ibn)
            for jk=npsp+1:size(att1,1)
                nbpoin=bndp(jk-npsp);
                rt=-(att1(jk,3)+att1(jk,4)) ...
                    *kapaFun(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +(att1(jk,1)*vecSp1(ppp(ipoin,1),ppp(ipoin,2),Tnow) ...
                    +att1(jk,2)*vecSp2(ppp(ipoin,1),ppp(ipoin,2),Tnow));
                acoe(ipoin,g2locNmlPbnd(nbpoin)+npoin)=acoe(ipoin,g2locNmlPbnd(nbpoin)+npoin) ...
                    +rt*thet*muFun(ppp(ipoin,1),ppp(ipoin,2),Tnow,dlt);
                               
            end             
          %  Fnum(ipoin)=Fnum(ipoin)+bndnmlVadd; 
        end
        
        if typPoints(ipoin)==2
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
            
            
            nor=ppp(ipoin,:); %if not unit normal direction, must change
     %      nBoundEq=nBoundEq+1;
            %zeroORnpoin
            if HRBFDQ==1
                % rt=0.0;
                for jk=1:n_pointPoint2(ipoin)    
                    acoe(npoin+g2locNmlPbnd(ipoin),pointsPoint2(ipoin,jk))=rd2(jk,1)*nor(1)+rd2(jk,2)*nor(2);
                end
                tmljm=tmpCell{ipoin,3};
                tmbptljm=tmpCell{ipoin,2};
                
                % xy1=ppp(ipoin,:);
                nd=n_pointPoint2(ipoin)+1;

                acoe(npoin+g2locNmlPbnd(ipoin),ipoin)=rd2(nd,1)*nor(1)+rd2(nd,2)*nor(2);
                
                npnb=size(pxynb,1);
            
                rt=0;
                for jk=1:npnb
%                     tm1=dfx1(pxynb(jk,1),pxynb(jk,2),Tnow)*pxynbnor(jk,1);
%                     tm2=dfy1(pxynb(jk,1),pxynb(jk,2),Tnow)*pxynbnor(jk,2);
%                     rt=rt+(rd2(jk+nd,1)*nor(1)+rd2(jk+nd,2)*nor(2))*(tm1+tm2);
                    acoe(npoin+g2locNmlPbnd(ipoin),npoin+g2locNmlPbnd(tmbptljm(jk)))=  ...
                        rd2(jk+nd,1)*nor(1)+rd2(jk+nd,2)*nor(2);

                end
                xy1=ppp(ipoin,:);
                tm1= dfx1(xy1(1),xy1(2),Tnow)*nor(1)+dfy1(xy1(1),xy1(2),Tnow)*nor(2);

                Fnum(npoin+g2locNmlPbnd(ipoin))=tm1;

                
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
%     diagacoe=zeros(npoin,npoin);
%    for jjmm=1:npoin
%        diagacoe(jjmm,jjmm)=1/acoe(jjmm,jjmm);
%    end
 %    diagacoe=pinv(diag(diag(acoe)));
% % 
% %   % unum(:,nStep+1)=pcg(acoe,Fnum,1e-6,100,diagacoe);
% %  % unum(:,nStep+1)=pinv(acoe)*Fnum;
% %  %   unum(:,nStep+1)=GuassFull(acoe,Fnum);
  %    unum(:,nStep+1)=gmres(acoe,Fnum,npoin,1e-12,150,diagacoe);

end
        
condAcoe=cond(acoe,2);        
        
uerr=unum(1:npoin,NtimeStep+1)-uexact(ppp(:,1),ppp(:,2),Tend);  
lrmserr=sqrt(sum(uerr.^2)/npoin);
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
plot3(ppp(:,1),ppp(:,2), unum(1:npoin,NtimeStep+1), 'b.','MarkerSize',20);
xlabel('x'); ylabel('y');
zlabel('u^h({\bf x}, T)');
grid on

figure(4)
plot3(ppp(:,1),ppp(:,2), abs(uerr), 'b.','MarkerSize',20);
zlabel('Error');
xlabel('x'); ylabel('y');
grid on

fid11=fopen('yuntu.plt','w');
nem=size(ttt,1);
fprintf(fid11, 'TITLE="u numerical solution"\n');
fprintf(fid11, 'VARIABLES="x","y","u","error"\n');
fprintf(fid11, 'ZONE N=%d,E=%d, F=FEPOINT, ET=TRIANGLE\n',npoin,nem);
for ij=1:npoin
    fprintf(fid11,'%f   %f   %f   %f\n',ppp(ij,1),ppp(ij,2),unum(ij,NtimeStep+1),uerr(ij));
end

for ij=1:nem
    fprintf(fid11,'%d  %d  %d\n', ttt(ij,1),ttt(ij,2),ttt(ij,3));
end

fclose(fid11);