%% q2
N=1000;
x=zeros(N,1);
y=zeros(N,1);
P=[0.8, 0.2;
    0.2, 0.8];
Pt=[0.8 0 0.2 0;
    0.8 0 0.2 0;
    0 0.2 0 0.2;
    0 0.2 0 0.2];
h1=1;
for i=1:N
    x(i)=-1+2*round(rand);
    if i~=1
        y(i)=x(i)+h1*x(i-1)+randn(1);
    end
end


%% q3
N=1000;
ite=10000;
M=3;
pi=[0.3,0.5,0.2];
lamda=[1,2,5];
% Observation:
y=zeros(N,1);
for n=1:N
    randpi=rand();
    if randpi<=0.2
        y(n)=exprnd(1);
    elseif randpi<=0.7
            y(n)=exprnd(1/2);
    else
        y(n)=exprnd(1/5);
    end
end
%EM algorithm
pi_m=[0.25,0.55,0.2];
lamda_m=[1.2, 2.5, 4.6];
w=zeros(N,M);
for t=1:ite
    for i=1:N
        %E step:
        w_denom=0;
        for jp=1:M
            w_denom=w_denom+pi_m(jp)*exppdf(y(i),1/lamda_m(jp));
        end
        for j=1:M
            w(i,j)=pi_m(j)*exppdf(y(i),1/lamda_m(j))/w_denom;
        end
    end
    %M step:
    w_sum=zeros(M,1);
    wy_sum=zeros(M,1);
    for j=1:M
        for i=1:N
            w_sum(j)=w_sum(j)+w(i,j);
            wy_sum(j)=wy_sum(j)+w(i,j)*y(i);
        end
        pi_m(j)=w_sum(j)/N;
        lamda_m(j)=w_sum(j)/wy_sum(j);
    end
end
pi_m
lamda_m