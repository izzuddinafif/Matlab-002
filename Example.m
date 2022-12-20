clc; clear; close;
format longG;
r=0.265;
x=350;
y=3.85;
N=30;
for k=1:N
    if k==1
        A(k)=0;
        B(k)=1;
        T(k)=(x*(((A(k)/x)^y-(log(r)/B(k)))^(1/y)))-A(k);
        T(k)=round(T(k), N);
        continue
    end
    g=k-1;
    ak(k)=g/(3*g+1);
    A(k)=A(k-1)+(T(k-1)*ak(k));
    bk(k)=(4*g+1)/(3*g+1);
    B(k)=B(k-1)*bk(k);
    T(k)=(x*(((A(k)/x)^y-(log(r)/B(k)))^(1/y)))-A(k);
    T(k)=round(T(k), N);
end
TT=zeros(1,N);
TT(1)=T(1);
for k=1:N
    if k>1
        TT(k)=T(k)+TT(k-1);
    end
end
L=1461; % 3 years
Cf=100;
a=100;
b=50;
kk=100;
t=0:0.0001:1499.9999;
hh=(y/x).*(t/x).^(y-1);
h0=@(t) (y/x).*(t/x).^(y-1);
D=integral(h0,0,L);
for i=1:N
    tk(i)=round(TT(i)*10000);
end
dd=ones(1,29);
for i=1:N
    if i==1
        d(i)=hh(tk(i+1))-hh(tk(i));
        dd(i)=d(i);
    elseif i>1 && i<30
        d(i)=hh(tk(i+1))-hh(tk(i));
        dd(i)=d(i)+dd(i-1);
        h(i-1)=hh(tk(i))-dd(i-1);
    end
end
for i=1:29
    if i==1
        s1(i)=dd(i)*(L-TT(i));
        s2(i)=a+b*dd(i);
        J(i)=Cf*abs(D-s1(i))+s2(i);
    elseif i>1
        s1(i)=dd(i)*(L-TT(i))+s1(i-1);
        s2(i)=a+b*dd(i)+s2(i-1);
        J(i)=Cf*abs(D-s1(i))+s2(i);
    end
end
plot(1:29,J)