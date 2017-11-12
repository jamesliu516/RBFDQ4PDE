
clear all
clf
global ppp ttt meshden
meshden=0.1;

generateCircle; % call mesh generation
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
%���������Ŀ

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
                    
                   
               
        







