generateCircle;

global p t;

nelem=size(t,1);
npoin=size(p,1);

pstart=zeros(npoin,1); 
pend=zeros(npoin,1);
threeP=3;
maxElem=20;
elemPoint=zeros(npoin,maxElem); %elements aroung a point

n_elemPoint=zeros(npoin,1); %% the number of element with common point
pointElem=[]; % point with elements
for iele=1:nelem
    ti=t(iele,:);
    for i=1:threeP
        if n_elemPoint(ti(i))==0            
            n_elemPoint(ti(i))=n_elemPoint(ti(i))+1;
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
                elemPoint(ti(i),n_elemPoint(ti(i)))=iele;
            end
        end
    end
end
                
                
                
                    
                    
            

        
    
clear ptmp
