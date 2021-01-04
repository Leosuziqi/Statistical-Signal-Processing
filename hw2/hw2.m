%% EECE562 Assignment 2
% Ziqi Su 21410161
%% q1
P=[0.2 0.3 0.1 0.4;
   0.15 0.35 0.05 0.45;
   0.5 0.25 0.12 0.13;
   0.23 0.27 0.13 0.37];  %Aperiodic adn irreducible
q=1/4;
c=zeros(4,1);
for m=1:4
    c(m)=max(P(m,:))/q;
end
Y=1; % assume it starts from state 1
N=10000;
realization=zeros(1,N);
for n=1:N
    U1=rand();
    X=floor(4*U1)+1;
    U=rand();
    while U>=P(Y,X)/(c(Y)*q)
        U1=rand();
        X=floor(4*U1)+1;
        U=rand();
    end
    Y=X;
    realization(n)=Y;
end
plot(1:N,realization);
xlim([0 100]) % for better visualization, I only show part of realization here
title("Silmulation of Markov Chain sequence");
%% q2
%%%%%%%% Estimating the transition probabilities %%%%%%%%
np=zeros(4,4);
for n=1:N-1
    switch realization(n)
    % Find state = one of the cases
        case 1
            switch realization(n+1)
                % Find next state = one of the cases
                case 1
                    np(1,1)= np(1,1)+1;
                case 2
                    np(1,2)= np(1,2)+1;
                case 3
                    np(1,3)= np(1,3)+1;
                case 4
                    np(1,4)= np(1,4)+1;
            end
        case 2
            switch realization(n+1)
                case 1
                    np(2,1)= np(2,1)+1;
                case 2
                    np(2,2)= np(2,2)+1;
                case 3
                    np(2,3)= np(2,3)+1;
                case 4
                    np(2,4)= np(2,4)+1;
            end
        case 3
            switch realization(n+1)
                case 1
                    np(3,1)= np(3,1)+1;
                case 2
                    np(3,2)= np(3,2)+1;
                case 3
                    np(3,3)= np(3,3)+1;
                case 4
                    np(3,4)= np(3,4)+1;
            end
        case 4
            switch realization(n+1)
                case 1
                    np(4,1)= np(4,1)+1;
                case 2
                    np(4,2)= np(4,2)+1;
                case 3
                    np(4,3)= np(4,3)+1;
                case 4
                    np(4,4)= np(4,4)+1;
            end
    end
end
estimated_p=zeros(4,4);
% Normalize estimated_p
for m=1:4
    for n=1:4
        estimated_p(m,n)=np(m,n)/sum(np(m,:));
    end
end

%%%%%%%% Estimating the stationary distribution of the Markov Chain %%%%%%%%
[V,D]=eig(estimated_p.');
distribution_est=zeros(4,1);
for m=1:4
     distribution_est(m)=V(m,1)/sum(V(:,1));
end
distribution_est

%%%%%%%% Estimating the stationary distribution from the transition probability matrix %%%%%%%%
[V,D]=eig(P.');
distribution_P=zeros(4,1);
for m=1:4
    distribution_P(m)=V(m,1)/sum(V(:,1));
end
% distribution_est is close to distribution_P, which proved that my realization and estimator is correct 
distribution_P

%% q3
mu=[0.3 0.2 0.1 0.35 0.05];
%reconstruct
mu1=sort(mu);
N=5;
P=zeros(N,N);
muk=zeros(N,N);
muk(1,:)=mu1;
beta=zeros(N,1);
alpha=ones(N,1);
% iteratively compute beta, alpha and mu 
for k=1:N-1
    beta(k+1)=1-muk(k,k);
    alpha(k+1)=1-muk(k,k)/beta(k+1);
    muk(k+1,:)=muk(k,:)/beta(k+1);
end
% compute P matrix
for m=1:N
    for n=1:N
        if m==n
            P(m,n)=0;
        end
        if m>n
            P(m,n)=prod(alpha(1:n))*(1-alpha(n+1));
        end
        if m<n
            P(m,n)=prod(alpha(1:m))*muk(m,n)/beta(n);
        end
    end
end

% Estimating the stationary distribution from the transition probability matrix 
[V,D]=eig(P.');
distribution_P=zeros(5,1);
for m=1:5
    distribution_P(m)=V(m,1)/sum(V(:,1));
end      
% Results should be close to the given mu=[0.3 0.2 0.1 0.35 0.05]   
sort(mu)
distribution_P'   
%% q4
% Problem 4 is written by hand and attached after problem 5.

%% q5
% (a)
p=[1/3 1/3 1/3];
N=10000;
x1=zeros(N,1);
for n=1:N
    i_star=randi(3);
    U1=rand();
    switch i_star
        case 1
            x1(n)=U1;
        case 2
            x1(n)=power(U1,1/3);
        case 3
            x1(n)=power(U1,1/5);
    end
end  
figure
hist(x1,50) 
title("Histogram of the empirical distribution function");
% compare with the actual distribution
figure
x_actual=linspace(0,1);
plot(x_actual,(1+3*x_actual.^2+5*x_actual.^4)/3)
title("Histogram of the actual distribution function");
% (b)
% Assume n=4, alpha=[1/4, 1/4, 1/4, 1/4];
p=[1/4, 1/4, 1/4, 1/4];
N=10000;
x2=zeros(N,1);
for n=1:N
    i_star=randi(4);
    U2=rand();
    switch i_star
        case 1
            x2(n)=U2;
        case 2
            x2(n)=power(U2,1/2);
        case 3
            x2(n)=power(U2,1/3);
        case 4
            x2(n)=power(U2,1/4);
    end
end  
figure 
hist(x2,50)   
title("Histogram of the empirical distribution function");
% compare with the actual distribution
figure
x_actual=linspace(0,1);
plot(x_actual,(1+2*x_actual+3*x_actual.^2+4*x_actual.^3)/4)
title("Histogram of the actual distribution function");
%% q6
% Problem 6 is attached at the end.

