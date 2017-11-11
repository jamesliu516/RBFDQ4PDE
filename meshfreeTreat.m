
clear all
clf
global p t meshden;
meshden=0.04;

generateCircle;

nelem=size(t,1);
npoin=size(p,1);

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
%包含点的数目

%pointElem=[]; % point with elements

for iele=1:nelem
    ti=t(iele,:);
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
        ti=t(elemPoint(ipoin,ielem));
        for i=1:threeP
            if ti(i) ~= ipoin
                if n_pointPoint(ipoin)==0
                    n_pointPoint(ipoin)=n_pointPoint(ipoin)+1;
                    pointsPoin(ipoin,n_pointPoint(ipoin))=ti(i);
                else
                    flag=0;
                    for j=1:n_pointPoint(ipoin)
                        if pointsPoin(ipoin,j)==ti(i)
                            flag=1;
                            break
                        end
                    end
                   
               
        







