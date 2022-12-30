clc; clear; close all;
format longG;
r=0.265;
x=350;
y=3.85;
N=7;
for k=1:N
    if k==1
        A(k)=0;
        B(k)=1;
        T(k)=(x*(((A(k)/x)^y-(log(r)/B(k)))^(1/y)))-A(k);
        continue
    end
    g=k-1;
    ak(k)=g/(3*g+1);
    A(k)=A(k-1)+(T(k-1)*ak(k));
    bk(k)=(4*g+1)/(3*g+1);
    B(k)=B(k-1)*bk(k);
    T(k)=(x*(((A(k)/x)^y-(log(r)/B(k)))^(1/y)))-A(k);
end
TT=zeros(1,N);
TT(1)=T(1);
for k=1:N
    if k>1
        TT(k)=T(k)+TT(k-1);
    end
end
L=1095; % days of lease
Cf=100;
a=100;
b=50;
t=0.0001:0.0001:L;
hh=(y/x).*(t/x).^(y-1);
h0=@(ts) (y/x).*(ts/x).^(y-1);
D=integral(h0,0,L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    tk(i)=round(TT(i)*10000);
end
dd=ones(1,N);
for i=1:N    
%{
 if i==1
        h(i)=B(i)*h0(A(i)+t(1));
        hk(i)=hh(tk(i));
        d(i)=hk(i)-h(i);
    elseif i>1
        h(i)=B(i)*h0(A(i)+t(1));
        hk(i)=hh(tk(i));
        d(i)=hk(i)-h(i);
 
%}
    if i==1
        h(i)=B(i)*h0(A(i)+t(1));
        p(i)=hh(tk(i));
        d(i)=p(i)-h(i);
    elseif i>1
        h(i)=B(i)*h0(A(i)+t(1));
        p(i)=hh(tk(i));
        d(i)=p(i)-h(i);
    end
end
for i=1:N
    if i==1
        s1(i)=d(i)*(L-TT(i));
        s2(i)=a+b*d(i);
        J(i)=Cf*abs(D-s1(i))+s2(i);
    elseif i>1
        s1(i)=d(i)*(L-TT(i))+s1(i-1);
        s2(i)=a+b*d(i)+s2(i-1);
        J(i)=Cf*abs(D-s1(i))+s2(i);
    end
end
plot(1:N,J) 
