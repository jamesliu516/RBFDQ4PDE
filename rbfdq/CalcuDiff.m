clear
syms x y nx ny c xm ym
f=((x-xm)*nx+(y-ym)*ny)/sqrt((x-xm)^2+(y-ym)^2+c);
f_x=diff(f,x);

f_y=diff(f,y);
f_x2=diff(f,x,2);
f_xy=diff(f_x,y);
f_y2=diff(f,y,2);


% f_x2 =3*(x - xm)^2*(nx*(x - xm) + ny*(y - ym))/(c + (x - xm)^2  ...
%   + (y - ym)^2)^(5/2) - (nx*2*(x - xm))/(c + (x - xm)^2 + (y - ym)^2)^(3/2) ...
%   - (nx*(x - xm) + ny*(y - ym))/(c + (x - xm)^2 + (y - ym)^2)^(3/2);


% f_x2 =3*(x - xm)^2*(nx*(x - xm) + ny*(y - ym))/(c + (x - xm)^2  ...
%   + (y - ym)^2)^(5/2) - (3*nx*(x - xm)+ny*(y - ym))  ...
%   /(c + (x - xm)^2 + (y - ym)^2)^(3/2);

%  f_y2xx = 3*(y - ym)^2*(nx*(x - xm) + ny*(y - ym))/(c + (x - xm)^2 + (y - ym)^2)^(5/2) - (nx*(x - xm) + 3*ny*(y - ym))/(c + (x - xm)^2 + (y - ym)^2)^(3/2)
%f_xy =(3*(2*x - 2*xm)*(2*y - 2*ym)*(nx*(x - xm) + ny*(y - ym)))/(4*(c + (x - xm)^2 + (y - ym)^2)^(5/2)) - (nx*(2*y - 2*ym))/(2*(c + (x - xm)^2 + (y - ym)^2)^(3/2)) - (ny*(2*x - 2*xm))/(2*(c + (x - xm)^2 + (y - ym)^2)^(3/2))
  
%f_xy11 =(3*(x - xm)*(y - ym)*(nx*(x - xm) + ny*(y - ym)))/(c + (x - xm)^2 + (y - ym)^2)^(5/2) - (nx*(y - ym))/(c + (x - xm)^2 + (y - ym)^2)^(3/2) - (ny*(x - xm))/(c + (x - xm)^2 + (y - ym)^2)^(3/2)
