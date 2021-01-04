%mf
TRANS_XM=[0.9833 0.0167;
    0.8632 0.1368];
EMIS_B=[0.9825 0.0175;
    0.17 0.83];
p=[0.995 0.005];
%mm
TRANS_XM=[0.9846 0.0154;
    0.8707 0.1293];
EMIS_B=[0.9825 0.0175;
    0.17 0.83];
p=[0.9803 0.0

T_aug= [0 p; zeros(size(TRANS_XM,1),1) TRANS_XM];
E_aug=[zeros(1,size(EMIS_B,2)); EMIS_B];
%[seq,states] = hmmgenerate(1000,T_aug,E_aug);

%likelystates = hmmviterbi(seq, T_aug, E_aug);

%sum(states==likelystates)/1000;
N=10000;
obs=zeros(1,N);
for j=1:N
    ran=rand();
    if ran>0.4
        obs(j)=1;
    else
        obs(j)=2;
    end
end

%[TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS,'maxiterations',30)

%PSTATES = hmmdecode(obs,T_aug,E_aug)
likelystates = hmmviterbi(obs, T_aug, E_aug);
%sum(PSTATES(2,:))/10000
%sum(PSTATES(3,:))/10000
AP_count=0;
MF_count=0;
for m=1:N
    if obs(m)==2
        AP_count=AP_count+1;
        if likelystates(m)==2
            MF_count=MF_count+1;
        end
    end
end
p29=MF_count/AP_count