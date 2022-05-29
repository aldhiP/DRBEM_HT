% DRBEM d^2u/dx^2 + d^2u/dy^2 + (x^2 + y^3)u = (6*x*y^2 + 2*y + 2*x^3) + (x^2 + y^3)*(x^3*y^2 +x^2*y)
% dengan syarat batas
% y=0 ======> u= 0
% x=1 ======> un= 3*y^2 + 2*y
% y=1 ======> un= 2*x^3 + x^2
% x=0 ======> u= 0
% dan Solusi Analitik: u= x^3 * y^2 + x^2 * y


clear;clc;
L=input('Masukkan jumlah titik interior (bilangan kuadratik)=');
s=sqrt(L);
x=zeros(1,L);
y=zeros(1,L);
for i=1:s
    for j=1:s
        x((i-1)*s+j)=j/(s+1);
        y((i-1)*s+j)=i/(s+1);
    end
end
u=[x;y];
fid=fopen('interiorHT.txt','wt');
fprintf(fid,'%8.6f  %8.6f \n',u);
fclose(fid);