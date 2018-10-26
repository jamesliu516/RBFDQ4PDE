
global ppp ttt  pointboun typPoints domain nelem
%pointboun: boundary node number
global n_pointPoint2 pointsPoint2
global n_pointPoint pointsPoint su2mesh filenmsu2
%meshden=0.025;
%domain=1;
switch domain
    case 1  
        if su2mesh==1
            filenmsu2='square11by11.su2';            
            %loadsu2mesh; 
             loadsu2mesh4ref;   
        else
            generateRectangle;
        end
    case 2
        if su2mesh ==1
         %     filenmsu2='circleu1sss.su2'; 
          filenmsu2='circleWithCircle.su2';
           %loadsu2mesh; 
            loadsu2mesh4ref;
        else
            generateCircle; % call mesh generation
        end
        
    case 33
        if su2mesh ==1
         %     filenmsu2='circleu1sss.su2'; 
          filenmsu2='part4star.su2';
           %loadsu2mesh; 
            loadsu2mesh4ref;  
        end
    otherwise
        warning('Unexpected demain type. No mesh created.');
end

%ttt: the element array����Ԫ������ʾ��Ԫ����Щ�����
%ppp: ������
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
%������ĵ�Ԫ��Ŀ
n_pointPoint=zeros(npoin,1); %% the number of point(support points) with common point
%���������Ŀ��һ�㵥Ԫ��
typPoints=zeros(npoin,1); %% the type of the node: 0 inner point, 1 Dirchlet Boundary point
                          %2 Neumann Boundary condition point
                          
typPoints(pointboun)=1; %1 boudary point  default Dirichlet boundary

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

%�ҵ���һ�㵥Ԫ�ĵ�
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

%�ҵ��ڶ��㵥Ԫ�ĵ�
n_pointPoint2=n_pointPoint;
pointsPoint2=pointsPoint;
onlyNearestNeighbor=0;
if onlyNearestNeighbor==0
for ipn=1:npoin   %��ipn
    for ip=1:n_pointPoint(ipn) %��ipn��Χ��һ������Ŀn_pointPoint(ipn)
        %n_elemPoint(pointsPoint(ipn,ip))��Χ��һ���ip�����м�����Ԫ���˵�
        for ielem=1:n_elemPoint(pointsPoint(ipn,ip))   
            %elemPoint(pointsPoint(ipn,ip),ielem)��Χ��һ���ip����������Ԫ�ĵ�ielem����Ԫ���
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
end

%clear n_pointPoint;
%clear pointsPoint;
clear elemPoint;

fprintf('+       Done mesh generation and neighbor point finding.      +\n');
                            
                            
                            
                            

        


                    
                   
               
        







