clc; clear; close all;
format longG;
x=1;
y=2;
m=0.5;
n=0.5;
GG=@(yy) 1-(1-exp(-(yy/n).^m));
tau=2;
G=integral(GG, tau, Inf);
L=5; %years of lease
Cf=100;
Cn=200;
Ct=300;
a=100;
b=50;
t=0.0001:0.0001:L;
hh=(y/x).*(t/x).^(y-1);
h0=@(ts) (y/x).*(ts/x).^(y-1);
D=integral(h0,0,L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=101;
