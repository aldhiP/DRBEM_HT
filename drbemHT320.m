% DRBEM d^2u/dx^2 + d^2u/dy^2 + (x^2 + y^3)u = (6*x*y^2 + 2*y + 2*x^3) + (x^2 + y^3)*(x^3*y^2 +x^2*y) 
% dengan syarat batas
% y=0 ======> u= 0
% x=1 ======> un= 3*y^2 + 2*y
% y=1 ======> un= 2*x^3 + x^2
% x=0 ======> u= 0
% dan Solusi Analitik: u= x^3 * y^2 + x^2 * y.



clear;clc;
fid=fopen('inputHT320.txt'); % 320 segmen
indat=fscanf(fid,'%g %g %g %g %g \n',[5,inf]); % 5 kolom
fclose(fid);
xb=indat(1,:);%nilai absis pangkal segmen (x(k))
yb=indat(2,:);%nilai ordinat pangkal segmen (x(k))
bt=indat(3,:);%boundary type, jika derivatif normal diketahui bt=1, jk tdk bt=0
bv1=indat(4,:);%nilai u pada batas, jika diket. Jika tidak, maka nilainya 0.
bv2=indat(5,:);%nilai du/dn pada batas, jika diket. Jika tidak, maka nilainya 0.
N=length(xb)-1;

l=zeros(1,N);
xmid=zeros(1,N);
ymid=zeros(1,N);
nx=zeros(1,N);
ny=zeros(1,N);

for i=1:N
    l(i)=sqrt((xb(i+1)-xb(i))^2+(yb(i+1)-yb(i))^2);
    xmid(i)=0.5*(xb(i+1)+xb(i));
    ymid(i)=0.5*(yb(i+1)+yb(i));
    nx(i)=(yb(i+1)-yb(i))/l(i);
    ny(i)=(xb(i)-xb(i+1))/l(i);
end


fid=fopen('interiorHT.txt'); % 400 titik interior
data=fscanf(fid,'%g %g \n',[2,inf]);
fclose(fid);
xint=data(1,:);
yint=data(2,:);
L=length(xint);
M=N+L; % 385 Total  titik

x=zeros(1,M);
y=zeros(1,M);
for i=1:M
    if i<=N
        x(i)=xmid(i);
        y(i)=ymid(i);
    else
        x(i)=xint(i-N);
        y(i)=yint(i-N);
    end
end

F1=zeros(M,N);
F2=zeros(M,N);
for n=1:M
    for k=1:N
        A=(l(k))^2;
        B=(-ny(k)*(xb(k)-x(n))+nx(k)*(yb(k)-y(n)))*2*l(k);
        E=(xb(k)-x(n))^2 + (yb(k)-y(n))^2;
        if n==k
            F2(n,k)=0;
            F1(n,k)=l(k)/(2*pi)*(log(l(k))+(1+B/(2*A))*log(1+B/(2*A))-B/(2*A)*log(B/(2*A))-1);
        else
            F2(n,k)=l(k)/(2*pi)*integral(@(t)f2(nx(k),xb(k),x(n),ny(k),yb(k),y(n),A,B,E,t),0,1);
            F1(n,k)=l(k)/(4*pi)*integral(@(t)f1(A,B,E,t),0,1);
        end
    end
end

% matrix \rho
rho=zeros(M);
for i=1:M
    for j=1:M
        rho(i,j)=1+(r(x(i),y(i),x(j),y(j)))^2+(r(x(i),y(i),x(j),y(j)))^3;
    end
end


psi=zeros(M);
for n=1:M
    for m=1:M
        if n<=N
            lambda=0.5;
        else
            lambda=1;
        end
        sum=0;
        for k=1:N
            sum=sum+F1(n,k)*(chiX(x(k),y(k),x(m),y(m))*nx(k)...
            +chiY(x(k),y(k),x(m),y(m))*ny(k))-F2(n,k)*chi(x(k),y(k),x(m),y(m));
        end
        psi(n,m)=lambda*chi(x(n),y(n),x(m),y(m))+sum;
    end
end
% Creating matriks \mu
mu=psi/rho; %psi*invers(rho)

% Processing
a=zeros(M);
b=zeros(M,1);
for n=1:M
    for k=1:N
        if n==k
            delta=1;
        else
            delta=0;
        end
        if bt(k)==0
            a(n,k)=-F1(n,k);
            b(n)=b(n)+(0.5*delta-F2(n,k)+mu(n,k)*f(x(k),y(k)))*bv1(k)-mu(n,k)*g(x(k),y(k));
        else
            a(n,k)=F2(n,k)-mu(n,k)*f(x(k),y(k))-0.5*delta;
            b(n)=b(n)+F1(n,k)*bv2(k)-mu(n,k)*g(x(k),y(k));
        end
    end
    for k=N+1:M
        if n==k
            delta=1;
        else
            delta=0;
        end
        a(n,k)=-mu(n,k)*f(x(k),y(k))-delta;
        b(n)=b(n)-mu(n,k)*g(x(k),y(k));
    end
end
z=a\b;



U=zeros(1,M);
Un=zeros(1,M);
An=zeros(1,M);
Er=zeros(1,M);

for i=1:N
    if bt(i)==0
        U(i)=bv1(i);
        Un(i)=z(i);
    else
        U(i)=z(i);
        Un(i)=bv2(i);
    end
end


for i=N+1:M
    U(i)=z(i);
    Un(i)=0;
end


for i=1:M
    An(i)=(x(i))^3 * (y(i))^2 + (x(i))^2 * y(i); %solusi analitik DRBEM di soal.
    Er(i)=abs(An(i)-real(U(i)));
end


t=[U;Un]; %U: Solusi numerik, Un: Solusi derivatif numerik
fid=fopen('U-dan-Un-HT_320.txt','wt');
fprintf(fid,'%8.6f  %8.6f \n',t);
fclose(fid);



u=[U;An;Er]; %U: Solusi numerik, An: Solusi Analitik, Er: Error
fid=fopen('Numerik_vs_Analitik_HT_320.txt','wt');
fprintf(fid,'%8.6f  %8.6f   %8.6f \n',u);
fclose(fid);