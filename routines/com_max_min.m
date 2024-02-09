function [A,T]=com_max_min(X_max,T_max,X_min,T_min)
%% Amplitude calculation in the PFF algorithm by combining maxima and minima
%
% input: 
% X_max: fitted local peaks
% T_max: corresponding time of the fitted local peaks
% X_min: fitted local valleys
% T_min: corresponding time of the fitted local valleys

% output:
% A: the amplitude (positive)
% T: corresponding time of the amplitude

%%
T=zeros(length(T_max)+length(T_min),1);
A=zeros(length(X_max)+length(X_min),1);
if length(T_max)<length(T_min) 
%     str1='min'; % minimum comes first
    for loop1=1:length(T_max)
        T(2*loop1-1)=T_min(loop1);
        T(2*loop1)=T_max(loop1);
        A(2*loop1-1)=-X_min(loop1);
        A(2*loop1)=X_max(loop1);
    end
    T(end)=T_min(end);
    A(end)=-X_min(end);
elseif length(T_max)>length(T_min)
%     str1='max'; % maximum comes first
    for loop1=1:length(T_min)
        T(2*loop1-1)=T_max(loop1);
        T(2*loop1)=T_min(loop1);
        A(2*loop1-1)=X_max(loop1);
        A(2*loop1)=-X_min(loop1);
    end
    T(end)=T_max(end);
    A(end)=X_max(end);
elseif length(T_max)==length(T_min)
    if T_max(1)<T_min(1)
%         str1='equal max';
        for loop1=1:length(T_min)
            T(2*loop1-1)=T_max(loop1);
            T(2*loop1)=T_min(loop1);
            A(2*loop1-1)=X_max(loop1);
            A(2*loop1)=-X_min(loop1);
        end
    else
%         str1='equal min';
        for loop1=1:length(T_max)
            T(2*loop1-1)=T_min(loop1);
            T(2*loop1)=T_max(loop1);
            A(2*loop1-1)=-X_min(loop1);
            A(2*loop1)=X_max(loop1);
        end
    end
end
 
end