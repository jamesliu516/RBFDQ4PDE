
global ppp ttt meshden pointboun
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
%meshden=0.025;

generateCircle; % call mesh generation
%ttt: the element array，单元数组显示单元有哪些点组成
%ppp: 点坐标
nelem=size(ttt,1);
npoin=size(ppp,1);

%pstart=zeros(npoin,1);
%pend=zeros(npoin,1);
threeP=3;
maxElems=20;
maxPoints=25;
elemPoint=zeros(npoin,maxElems); %elements aroung a point
pointsPoint=zeros(npoin,maxPoints);

n_elemPoint=zeros(npoin,1); %% the number of element with common point
%包含点的单元数目
n_pointPoint=zeros(npoin,1); %% the number of point(support points) with common point
%包含点的数目第一层单元的
typPoints=zeros(npoin,1); %% the type of the node: 0 inner point, 1 Dirchlet Boundary point
                          %2 Neumann Boundary condition point
                          
typPoints(pointboun)=1;

%pointElem=[]; % point with elements

for iele=1:nelem
    ti=ttt(iele,:);
    for i=1:threeP
        if n_elemPoint(ti(i))==0
            n_elemPoint(ti(i))=n_elemPoint(ti(i))+1;
            if n_elemPoint(ti(i)) > maxElems
                disp('the number of element with common point more than maxElem');
                return;
            end
            elemPoint(ti(i),n_elemPoint(ti(i)))=iele;
        else
            flag=0;
            for j=1:n_elemPoint(ti(i))
                if elemPoint(ti(i),j)==iele
                    flag=1;
                    break
                end
            end
            if flag==0
                n_elemPoint(ti(i))=n_elemPoint(ti(i))+1;
                if n_elemPoint(ti(i)) > maxElems
                    disp('the number of element with common point more than maxElem');
                    return;
                end
                elemPoint(ti(i),n_elemPoint(ti(i)))=iele;
            end
        end
    end
end

%找到第一层单元的点
for ipoin=1:npoin
    for ielem=1:n_elemPoint(ipoin)
        ti=ttt(elemPoint(ipoin,ielem),:);
        for i=1:threeP
            if ti(i) ~= ipoin
                if n_pointPoint(ipoin)==0
                    n_pointPoint(ipoin)=n_pointPoint(ipoin)+1;
                    pointsPoint(ipoin,n_pointPoint(ipoin))=ti(i);
                else
                    flag=0;
                    for j=1:n_pointPoint(ipoin)
                        if pointsPoint(ipoin,j)==ti(i)
                            flag=1;
                            break
                        end
                    end
                    
                    if flag==0                        
                        n_pointPoint(ipoin)=n_pointPoint(ipoin)+1;
                        if n_pointPoint(ipoin)>maxPoints                          
                            disp('the number of points with common point more than maxPoints');
                            return;
                        end                        
                        pointsPoint(ipoin,n_pointPoint(ipoin))=ti(i);
                    end
                end
            end
        end
    end
end

%找到第二层单元的点
n_pointPoint2=n_pointPoint;
pointsPoint2=pointsPoint;
for ipn=1:npoin   %点ipn
    for ip=1:n_pointPoint(ipn) %点ipn周围第一层点的数目n_pointPoint(ipn)
        %n_elemPoint(pointsPoint(ipn,ip))周围第一层第ip个点有几个单元共此点
        for ielem=1:n_elemPoint(pointsPoint(ipn,ip))   
            %elemPoint(pointsPoint(ipn,ip),ielem)周围第一层第ip个点所共单元的第ielem个单元序号
            ti=ttt(elemPoint(pointsPoint(ipn,ip),ielem),:);
            flag=0;
            
            for ie=1:n_elemPoint(ipn)
                if elemPoint(pointsPoint(ipn,ip),ielem)==elemPoint(ipn,ie)
                    flag=1;
                    break;
                end
            end
            
            if flag==0
                for itm0=1:threeP
                    flag1=0;
                    for itm=1:n_pointPoint2(ipn)                       
                        if pointsPoint2(ipn,itm) == ti(itm0)
                            flag1=1;
                            break;
                        end
                    end
                    
                    
                    if flag1==0
                        n_pointPoint2(ipn)=n_pointPoint2(ipn)+1;
                        if n_pointPoint2(ipoin)>maxPoints
                            disp('the number of points2 with common point more than maxPoints');
                            return;
                        end                       
                        pointsPoint2(ipn,n_pointPoint2(ipn))=ti(itm0);
                    end
                end
            end
        end
    end
end

clear n_pointPoint;
clear pointsPoint;
clear elemPoint;

                            
                            
                            
                            

        


                    
                   
               
        







