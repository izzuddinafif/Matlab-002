%set parameter;
Cf = 10; 
A0 = 8;
fc = 10;
vc = 4;
L = 15;

%set parameter stage 1;
k = 10; %jumlah PM action -> indexing 
lb = 0; %batas bawah dari t
ub = L-(vc/Cf); %batas atas dari t
t = lb + (ones(1,k)).*(rand(1,k)).*(ub-lb); %t itu terletak antara 0 sampai L-b/Cf
lambda_0 = [0.999]; %buat initial lambda
lambda_t = [0.998 0.986 0.976 0.965 0.761 0.644 0.632 0.521 0.420 0.420]; %nilai lambda

%ngitung lambda t min
for j=1:k
    if j==1
        lambda_tmin(j) = lambda_t(j)-lambda_0;
    else
        lambda_tmin(j) = lambda_t(j)-lambda_t(j-1);
    end
end

do = optimvar('do',1,k,'LowerBound',0); %set variabel yang akan dioptimasi
%fungsi objektif
for j = 1:k
    obj = Cf.*A0 + fc.*k - Cf.*(sum(do(j).*(L-t(j)-(vc./Cf))));
end
prob1 = optimproblem('Objective',obj,'ObjectiveSense','min'); %membentuk rangkaian problem
for j = 1:k
    cons = do(j) <= lambda_tmin(j);
end
prob1.Constraints.cons = cons; %memasukkan constrain ke rangkaian problem

do_star = cell2mat(struct2cell(solve(prob1))) %solve probelm, tapi format jawabannya diubah dalam bentuk matrix


%========================STAGE 2;
L_bar = L-(vc/Cf); %mencari l bar dr rumus
t2 = optimvar('t2',1,k,'LowerBound',0); %tj* sebagai variabel keputusan
%fungsi objektif
for j = 1:k
    if j == 1
        obj2 = Cf.*(A0-(L-(vc./Cf).*lambda_t(k)+(sum(t2(j).*(lambda_t(j)-lambda_0)))))+k*fc;
    else
        obj2 = Cf.*(A0-(L-(vc./Cf).*lambda_t(k)+(sum(t2(j).*(lambda_t(j)-lambda_t(j-1))))))+k*fc;
    end
end
prob2 = optimproblem('Objective',obj2,'ObjectiveSense','min'); %membuat rangkaian problem
%constrain
for j = 1:k
    cons = t2(j) <= L_bar;
end
prob2.Constraints.cons = cons; %constrain dimasukkan ke rangkaian problem

t_star = cell2mat(struct2cell(solve(prob2))) %solve problem dengan format matrix


%========================STAGE 3;
%memasukkan do* dan tj* ke pers J cost, 
%lalu dicari yang paking minimum di index ke brp untuk jadi k*
for j = 1:k
    J_cost(j) = Cf.*(A0-(sum(do_star(j).*(L-t_star(j)))))+(sum(fc+vc.*do_star(j)));
end

J_cost
