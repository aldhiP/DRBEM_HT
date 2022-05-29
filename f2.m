function y=f2(nx,xb,xmid,ny,yb,ymid,A,B,E,t)
y=((xb-xmid)*nx+(yb-ymid)*ny)./(A*t.^2+B*t+E);