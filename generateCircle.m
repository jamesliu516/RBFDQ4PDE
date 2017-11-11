%Generate circl example.

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
clear 
clc
global p t meshden
rand('state',1); % Always the same results
set(gcf,'rend','z');

fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));

fprintf('Uniform Mesh on Unit Circle\n');
echo on
ra=1;
fd=@(p) sqrt(sum(p.^2,2))-ra;
[p,t]=distmesh2d(fd,@huniform,meshden,[-ra,-ra;ra,ra],[]);
echo off
fstats(p,t);
fprintf('(Done mesh generation.)\n\n')

