Ra=2;
rr=1;
x(1)=0;
y(1)=Ra;

x(2)=cos(126*pi/180)* rr;
y(2)=sin(126*pi/180)* rr;

x(3)=-cos(18*pi/180)*Ra;
y(3)=sin(18*pi/180)*Ra;

x(4)=cos((126+72)*pi/180)*rr;
y(4)=sin((126+72)*pi/180)*rr;

x(5)=cos(234*pi/180)*Ra;
y(5)=sin(234*pi/180)*Ra;

x(6)=cos((126+72+72)*pi/180)*rr;
y(6)=sin((126+72+72)*pi/180)*rr;

x(7)=cos((234+72)*pi/180)*Ra;
y(7)=sin((234+72)*pi/180)*Ra;

plot(x,y)
