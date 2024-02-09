function [Freq3,damp3,RespDecom,RespPred]=ERA_SISO_decay(data_0,nSRM,ConstOrder,fs,PointPredict)
%% System Reconstruction and Prediction based on Eigensystem Realization Algorithm (ERA)
%
% input:
% data_0: the data chosen for identification (must be row vector)
% nSRM: the order of the reconstructed system
% ConstOrder: the orders that need to be constructed, e.g. [1 2 4]
% fs: sampling frequency
% PointPredict: integer, the data points that are used for system reconstruction
%
% output:
% Freq3: the identified frequencies based on the chosen data (usually not needed)
% damp3: the identified damping based on the chosen data (usually not needed)
% RespDecom: the decomposed data of the chosen data (usually not needed)
% RespPred: the data points that are predicted
%
%    Wei Chen
%    meshiawei@tongji.edu.cn 
%    Version: 20200521

%%
RevFactor=1; % set revise factor if needed for experimental data
nRow=size(data_0,1);
Hr=2*nSRM+2;    % the row of Hankel matrix; which can be modified if necessary.
Hc=nRow-Hr;    %the column of Hankel matrix;
% Construct the  Hankel matrix 0
H_0=zeros(Hr,Hc);
for loop1=1:Hr
    for loop2=1:Hc
        H_0(loop1,loop2)=data_0(loop1+loop2-1,:);
    end
end

% Construct the  Hankel matrix 1
H_1=zeros(Hr,Hc);
for loop3=1:Hr
    for loop4=1:Hc
        H_1(loop3,loop4)=data_0(loop3+loop4,:);
    end
end

% conduct svd for the H_0
[U_svd1,S_svd1,V_svd1]=svd(H_0);
U_svd2=U_svd1(:,1:2*nSRM);
S_svd2=S_svd1(1:2*nSRM,1:2*nSRM);
V_svd2=V_svd1(:,1:2*nSRM);
% system matric, control matric, and the observe matric
As=S_svd2^(-1/2)*U_svd2.'* H_1*V_svd2*S_svd2^(-1/2);
[EigVec,EigVal]=eig(As);
EO=[1,zeros(1,Hr-1)];
C_SRM=EO*U_svd2*S_svd2^(1/2);
EC=[1,zeros(1,Hc-1)].';
B_SRM=S_svd2^(1/2)*(V_svd2.')*EC; % initial state variable 

% modal parameters
EigValrow=[diag(EigVal)].';
lmd1=log(EigValrow)*fs;    % calculate the eigenvalues£¬lmt=log(eigenvalues)/dt
Freq1=abs(lmd1)/(2*pi);       % natural frequency£¬Freq=abs(lmt)
damp1=sqrt(1./(((imag(log(EigValrow))./real(log(EigValrow))).^2)+1));    %calculate the damping ratio

% sequence the frequency
AllInform=[Freq1;damp1;lmd1;EigVec];
for loop5=1:size(Freq1,2)-1
    for loop6=loop5+1:size(Freq1,2)
        if AllInform(1,loop5)<=AllInform(1,loop6)
        else
           tempVec=AllInform(:,loop5); AllInform(:,loop5)=AllInform(:,loop6);AllInform(:,loop6)=tempVec;
        end
    end
end
EigVec2=AllInform(4:end,:);Freq2=AllInform(1,:);damp2=AllInform(2,:);lmd2=AllInform(3,:);

%% Response Seperation
ModalRespon=inv(EigVec)*S_svd2^(1/2)*V_svd2.';
ncount=1;
for loop7=1:size(Freq2,2)-1
    if Freq2(loop7)==Freq2(loop7+1)
        Freq3(ncount,:)=Freq2(loop7);
        damp3(ncount,:)=damp2(loop7);
        Stateresponse=real(EigVec(:,[loop7,loop7+1])*ModalRespon([loop7,loop7+1],:));
        RespDecom(ncount,:)=C_SRM*Stateresponse;
        ncount=ncount+1;
    end
end

%% Response Prediction
eigenspace=EigVec2\B_SRM;
ncount=1;
for loop9=ConstOrder
lmt3(2*ncount-1,2*ncount-1)=exp(lmd2(2*loop9-1)/fs);   % the
lmt3(2*ncount,2*ncount)=exp(lmd2(2*loop9)/fs);   % the
EigVec3(:,2*ncount-1)=EigVec2(:,2*loop9-1);
EigVec3(:,2*ncount)=EigVec2(:,2*loop9);
TarEig(2*ncount-1,:)=eigenspace(2*loop9-1);
TarEig(2*ncount,:)=eigenspace(2*loop9);
ncount=ncount+1;
end
for loop8=1:PointPredict
    if loop8==1
    lmtuse=lmt3;
    lmt4=diag((diag(lmt3)).^(-1));
    else
    lmtrev=log(lmtuse)*fs;    % calculate the eigenvalues£¬lmt=log(eigenvalues)/dt
    lmtrev=real(lmtrev)*RevFactor+1j*imag(lmtrev);
    lmtuse=exp(lmtrev/fs);
    lmt4=lmt4*diag((diag(lmtuse)).^(-1));
    end
    StatePred(:,loop8)= real(EigVec3*lmt4*TarEig);
end

RespPred=[C_SRM*StatePred].';
% sr=RespPred
end

