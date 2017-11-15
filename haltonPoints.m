%halton ponts in 
clf
rh=1;
nump=256*4;
ah=haltonset(2,'Skip',1e3,'Leap',1e2);
ah1=net(ah,nump);

plot(ah1(:,1),ah1(:,2),'o')
