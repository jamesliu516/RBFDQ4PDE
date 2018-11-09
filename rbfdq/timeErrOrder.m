%time error order

timestep=[1.0/50; 1.0/100; 1.0/200; 1.0/400];
errl2=[5.7446e-4; 2.4377e-4; 9.7826e-5; 3.3469e-5];

a=ones(length(timestep),2);
a(:,2)=-log(timestep);
b=-log(errl2);

jie=a\b;

%plot(-log(timestep), b,':o', -log(timestep), 3.2-log(timestep),'b-', ...
 %   -log(timestep), -2*log(timestep),'r-')
% include second layer grid
c2=[2.5;5; 8; 10; 12; 15; 18; 20; 22;25;30;35;40];
c2l2err=[4.1089e-4;3.5034e-5;5.6596e-5; 7.6915e-5; 8.8299e-5;...
    9.7826e-5;1.0308e-4;1.0541e-4; 1.0713e-4; 1.0898e-4; 1.1096e-4;1.1226e-4; 1.1301e-4];

%hold on
c2L=[2.5;5; 8; 10; 12; 15; 18; 20; 22;25;30;35;38;40;42;45;50;55;65;75;95];
c2l2errL=[2.2924e-2; 1.2498e-2; 8.1159e-3; 6.5891e-3; 5.5517e-3; 4.4982e-3; ...
    3.7866e-3; 3.4281e-3; 3.1333e-3;  2.7779e-3; 2.3409e-3;...
    2.0271e-3; 1.8780e-3;1.7908e-3; 1.7119e-3;1.6064e-003;1.4586e-003;...
    1.3374e-003;1.1506e-003; 1.0132e-003;8.2480e-04];

numbb=length(c2);
plot(c2L(1:numbb),log10(c2l2errL(1:numbb)),'-o',c2,log10(c2l2err),'-+','LineWidth',2, 'MarkerSize',7)
xlabel('c^2')
ylabel('log_{10}(L^2 error)')

legend('One layer','Two layer')
min1=10000;
max1=0;
for ij=1:size(n_pointPoint2,1)
    if ij~=pointboun
        if n_pointPoint2(ij) < min1
            min1=n_pointPoint2(ij);
        end
        
        if n_pointPoint2(ij) > max1
            max1=n_pointPoint2(ij);
        end 
    end
end
min1
max1
        