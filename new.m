% Numerical Example
datetime
clear; close all; clc;
format longG
N=10;
x=350;
y=4;
dcp=2;
drp=[10,50,100,500];
Rth=[0.313,0.289,0.277,0.265];
Tc=4;
Tp=2;
r=0.001:0.001:1;
h=@(zs) (y/x)*(zs/x).^(y-1);
no=1; % !!!!!!!!! ganti ini untuk pilih drp (1-4; 1=10, 2=50, dst)
% Step 1: Compute solutions Tk of Equation (11). Theses solutions are functions of reliability threshold Rth
Rrth=0.31;
for k=1:N
    g=k-1;
    if k==1
        A(k)=0;
        B(k)=1;
        T(k)=(x*(((A(k)/x)^y-(log(Rrth)/B(k)))^(1/y)))-A(k);
        continue
    end
    a(k)=g/(3*g+2);
    A(k)=A(k-1)+(T(k-1)*a(k));
    b(k)=(2*g+3)/(g+2);
    B(k)=B(k-1)*b(k);
    T(k)=(x*(((A(k)/x)^y-(log(Rrth)/B(k)))^(1/y)))-A(k);
end
% Step 2: Use solutions in Step 1 to solve Equation (23) with respect to Rth
drp=drp(no);
for i=1:N
    for j=1:N
        if j==1
            s(j)=1/(B(j)*0.00001*(A(j)+T(j)));
            continue
        end
        s(j)=s(j-1)+(1/(B(j)*0.00001*(A(j)+T(j))));
    end  
    ava(i)=s(i)/((i-1)*dcp+drp+s(i));
end
% Step 3: Based on Steps 1 & 2, choose N to maximise the function given by Equation (24). From the plot, N=7 is the optimal value of N.
 % cari av maks dan N optimal
[M,I] = max(ava);
disp("Availability = %");
disp(M*100);
disp("Number of maintenance = ");
disp(I);
for i=1:1000
    for k=1:N
        if k==1
            A(k)=0;
            B(k)=1;
            T(k)=(x*(((A(k)/x)^y-(log(r(i))/B(k)))^(1/y)))-A(k);
            R=@(t) exp(-B(k)*integral(h,A(k),t+A(k),'ArrayValued',true));
            ss(k)=integral(R,0,T(k),'ArrayValued',true);
            %av(i,k)=ss(k)/((k-1)*dcp+drp+ss(k));
            av(i,k)=ss(k)/(((k-1)*(Tc*(1-r(i))+Tp*r(i)))+drp+ss(k));
            continue
        end
        g=k-1;
        ak(k)=g/(3*g+2);
        A(k)=A(k-1)+(T(k-1)*ak(k));
        bk(k)=(2*g+3)/(g+2);
        B(k)=B(k-1)*bk(k);
        T(k)=(x*(((A(k)/x)^y-(log(r(i))/B(k)))^(1/y)))-A(k);
        R=@(t) exp(-B(k)*integral(h,A(k),t+A(k),'ArrayValued',true));
        ss(k)=ss(k-1)+integral(R,0,T(k),'ArrayValued',true);
        %av(i,k)=ss(k)/((k-1)*dcp+drp+ss(k));
        av(i,k)=ss(k)/(((k-1)*(Tc*(1-r(i))+Tp*r(i)))+drp+ss(k));
    end
[Av(i),n(i)]=max(av(i,:));
end
[AA,V]=max(Av);
rr=V/100;
disp("Availability = %");
disp(AA*100);
disp("Rth = %")
disp(rr*100);
grid on
% Step 4: From results of Steps 1 & 2 together with the optimal number N from Step 3, compute Tk (k = 1, . . . , N),
clear T A a B b
Rth=Rth(no);
for k=1:I
    g=k-1;
    if k==1
        A(k)=0;
        B(k)=1;
        T(k)=(x*(((A(k)/x)^y-(log(rr)/B(k)))^(1/y)))-A(k);
        continue
    end
    a(k)=g/(3*g+2);
    A(k)=A(k-1)+(T(k-1)*a(k));
    b(k)=(2*g+3)/(g+2);
    B(k)=B(k-1)*b(k);
    T(k)=(x*(((A(k)/x)^y-(log(rr)/B(k)))^(1/y)))-A(k);
end
disp("T1, ... , Tk = ")
disp(T(1:I))
%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=I;
for k=1:N
    if k==1
        A(k)=0;
        B(k)=1;
        T(k)=(x*(((A(k)/x)^y-(log(Rth)/B(k)))^(1/y)))-A(k);
        continue
    end
    g=k-1;
    ak(k)=g/(3*g+1);
    A(k)=A(k-1)+(T(k-1)*ak(k));
    bk(k)=(4*g+1)/(3*g+1);
    B(k)=B(k-1)*bk(k);
    T(k)=(x*(((A(k)/x)^y-(log(Rth)/B(k)))^(1/y)))-A(k);
end
TT=zeros(1,N);
TT(1)=T(1);
for k=1:N
    if k>1
        TT(k)=T(k)+TT(k-1);
    end
end
L=1460; %days of lease
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
for i=1:N 
    if i==1
        h1(i)=B(i)*h0(A(i)+t(1));
        p(i)=hh(tk(i));
        d(i)=p(i)-h1(i);
    elseif i>1
        h1(i)=B(i)*h0(A(i)+t(1));
        p(i)=hh(tk(i));
        d(i)=p(i)-h1(i);
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
figure(1)
plot(1:N,J) 
grid on
figure(2)
plot(1:N,ava(1:N))
grid on
datetime