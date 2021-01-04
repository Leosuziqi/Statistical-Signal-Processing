%% q2
N=1000;
y=zeros(N,1);
y(1)=1;
y(2)=2;
a1=-1;
a2=-0.1;
for k=3:N
    wk=normrnd(0,1);
    y(k)=a1*y(k-1)+a2*y(k-2)+wk;
end
plot(1:N,y)

%% q3
% Asymptotically stationary
for N=1000:1000:3000
    %simulate AR model as q2
    y=zeros(N,1);
    y(1)=1;
    y(2)=2;
    a1=-1;
    a2=-0.1;
    for k=3:N
        wk=normrnd(0,1);
        y(k)=a1*y(k-1)+a2*y(k-2)+wk;
    end
    % LS estimator
    pha=zeros(N-2,2);
    Yout=zeros(N-2,1);
    Yout=y(3:N);
    pha(:,1)=y(2:N-1);
    pha(:,2)=y(1:N-2);
    theta=inv(pha'*pha)*pha'*Yout;
    %mean square error  
    accuracy(N/1000)=immse(theta,[a1;a2]);
end
figure
scatter(1000:1000:3000,accuracy)
title("Accuracy of the LS  estimate with increasing data length");

% Non-asymptotically stationary
%When the model is not asymptotically stationary, y(N) will go to
%infinite.Infinite values cannot be calculated here.

%% q4
%(a)
N=1000;
theta(1)=1;
eplist=[0.1,0.01,0.001];
for ep=eplist
    for k=1:N
        pha=normrnd(0,1);
        v=normrnd(0,1);
        wk=normrnd(0,1);
        theta(k+1)=theta(k)+ep*wk;
    end
    plot(1:N,theta(1:N));
    hold on
end
legend('ep=0.1','ep=0.01','ep=0.001')

%% (b)
N=1000;
pha=zeros(N,1);
theta=zeros(N,1);
v=zeros(N,1);
y=zeros(N,1);
theta(1)=0.6;
ep=0.0001;
for k=1:N
        pha(k)=normrnd(0,1);
        v(k)=normrnd(0,1);
        wk(k)=normrnd(0,1);
        theta(k+1)=theta(k)+ep*wk(k);
        y(k)=pha(k)*theta(k)+v(k);
end
% RLS
theta_est(1)=1;
rou=0.999;
p=zeros(N,1);
p(1)=0.1;
for k=1:N-1
    alpha=rou^(N-k-1);
    theta_est(k+1)=theta_est(k)+(p(k)*pha(k+1))/(1/alpha+pha(k+1)'*p(k)*pha(k+1))*(y(k+1)-pha(k+1)'*theta_est(k));
    p(k+1)=p(k)-(p(k)*pha(k+1)*pha(k+1)'*p(k))/(1/alpha+pha(k+1)'*p(k)*pha(k+1));
end
theta(N)
theta_est(N)




