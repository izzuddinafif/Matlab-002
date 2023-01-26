clc; clear all; close all;
a=1;
b=2;
L=5;
Cf=100;
A=100;
B=50;
LL=L-(B/Cf);
t=0.1:0.1:L;
h=@(x) (b/a)*(x/a).^(b-1);
hh=(b/a).*(t/a).^(b-1);
D=integral(h,0,L);
for j=1:10
    if j==1
        T(j)=0;
        v(j)=0;
        continue;
    end
    v(j)=(b-1)/(b-v(j-1)^(b-1));
    T(j)=v(j)*LL;
end
j=1;
for d=0.01:0.01:1
    if d==0.01
        s1(j)=d*(L-T(9));
        s2(j)=a+b*d;
        J(j)=Cf*abs(D-s1(j))+s2(j);
        j=j+1;
        continue;
    end
    s1(j)=d*(L-T(9))+s1(j-1);
    s2(j)=a+b*d+s2(j-1);
    J(j)=Cf*abs(D-s1(j))+s2(j);
    j=j+1;
end