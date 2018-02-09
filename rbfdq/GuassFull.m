function  x = GuassFull( A,B )
%说明A是方阵，b是列向量
%example:
%A=[1,2,3;4,8,9;10,11,12];
%B = [2,3,5];
%x =GuassFull(A,B');
% Created By wangxiaole 14.9.18
% This function use Guass methods to calculate Linear Equations
% Choose the maximum element in all Matrix A;
[m,n] = size(A);
if rank(A) == m
    if  m == 1 %只有一个方程
        if A ~=0
            x = B/A;
        else
            x = inf;
        end
    else
        
        MaxEle = max(max(abs(A)));
        [m1,n1] = find(abs(A) == MaxEle);
        Row = A(:,n1(1));
        A(:,n1(1))=[];
        A = [Row,A];
        Column  =A(m1(1),:);
        A(m1(1),:) =[];
        A = [Column;A];
        
        B_Column = B(m1(1));
        B(m1(1)) = [];
        B =[B_Column ;B];
        
        for i = 2:m
            temp  = A(i,1);
            
            A(i,:) = A(i,:) + -1*temp/A(1,1)*A(1,:);
            B(i)= -1*temp/A(1,1)*B(1) +B(i);
        end
        A1 = A;
        B1 = B;
        
        A1(1,:)=[];
        A1(:,1)=[];
        B1(1)= [];
        x = GuassFull(A1,B1);
        tempx = (sum(-1*A(1,2:end).*x') +B(1))/A(1,1);
        if(x == inf)
            x = inf;
        else
            if Row == 1
                x = [tempx;x];
            else
                x = [x(1:n1(1)-1);tempx;x(n1(1):end)];
            end
        end
    end
else
    x = inf;
end

end