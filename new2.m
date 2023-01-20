clc; clear all; close all;
a=1;
b=2;
L=5;
Cf=100;
A=100;
B=50;
LL=L-(b/Cf);
t=0.1:0.1:L;
h=@(x) (b/a)*(x/a).^(b-1);
hh=(b/a).*(t/a).^(b-1);
D=integral(h,0,L);
for j=1:101
    if j==1
        t(j)=0;
        v(j)=0;
        continue;
    end
    v(j)=(b-1)/(b-v(j-1)^(b-1));
    t(j)=v(j)*LL;
end