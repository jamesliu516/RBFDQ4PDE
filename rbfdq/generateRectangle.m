

%Generate Rectangle example.

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
cd ../
global ppp ttt meshden pointboun racLow racHigh
%meshden=0.05; %defaut

pointboun=[];
rand('state',1); % Always the same results
set(gcf,'rend','z');

fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));


fprintf('Square, with size function point and line sources\n');

%[-1,1]X[-1,1]
fd=@(p) drectangle(p,racLow(1),racHigh(1),racLow(2),racHigh(2));

%[0,1]X[0,1]
%fd=@(p) drectangle(p,0,1,0,1);
% fh=@(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
%                 0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
%[-1,1]X[-1,1]
[ppp,ttt]=distmesh2d(fd,@huniform,meshden, ...
    [racLow(1),racLow(2);racHigh(1),racHigh(2)], ...
    [racLow(1),racLow(2);racHigh(1),racLow(2);racLow(1),racHigh(2); ...
    racHigh(1),racHigh(2)]);
%[0,1]X[0,1]
%[ppp,ttt]=distmesh2d(fd,@huniform,meshden,[0,0;1,1],[0,0;1,0;0,1;1,1]);
fstats(ppp,ttt);
%fprintf('(press any key)\n\n'); pause                  
                  
% fprintf('Uniform Mesh on Unit Circle\n');
% echo on
% ra=1;
% fd=@(p) sqrt(sum(p.^2,2))-ra;
% [ppp,ttt]=distmesh2d(fd,@huniform,meshden,[-ra,-ra;ra,ra],[]);
% echo off
% fstats(ppp,ttt);

eboun=boundedges(ppp,ttt);
neb=size(eboun,1);
for i=1:neb
    pointboun=[pointboun eboun(i,:)];
end
pointboun=sort(pointboun);
pointboun=pointboun';
pointboun=unique(pointboun);
clear eboun;
fprintf('(Done mesh generation.)\n\n')
cd rbfdq
