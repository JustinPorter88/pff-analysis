function [AMP_avg,FRE_avg,FRE_upper,FRE_lower,DAM_avg,DAM_upper,DAM_lower]=AvgCiC(FRE,AMP,DAM)
%% Calculation of averaged frequency and damping ratio with 95% confidence intervals
%
% input: (with the same length)
% FRE: frequency (column vector)
% AMP: amplitude (column vector)
% DAM: damping ratio (column vector)
%
% output:
% AMP_avg: averaged amplitude
% FRE_avg: averaged frequency
% FRE_upper: upper boundary of frequency
% FRE_lower: lower boundary of frequency
% DAM_avg: averaged damping ratio
% DAM_upper: upper boundary of damping ratio
% DAM_lower: lower boundary of damping ratio

%% determine the intersection of the amplitude
A_max = max(AMP{1});
for j = 2:length(AMP)
  A_max = min(A_max,max(AMP{j})) ;
end
A_min = min(AMP{1});
for j = 1:length(AMP)
  A_min = max(A_min,min(AMP{j})) ;
end
AMP_all_0 = [] ;
for j = 1:length(AMP)
  idx1 = find(AMP{j}<=A_max,1,'first') ;
  idx2 = find(AMP{j}>=A_min,1,'last') ;
  AMP_all_0 = [AMP_all_0; AMP{j}(idx1:idx2)]; %#ok<*AGROW>
end

%% sort the amplitude and return only unique values
AMP_all=unique(AMP_all_0);
clear AMP_all_0 A_min A_max idx1 idx2

%% interpolation for FRE and DAM according to AMP_all
FRE_int = zeros(length(AMP_all),length(AMP)) ;
DAM_int = FRE_int ;
for j=1:length(AMP)
    FRE_int(:,j)=interp1(AMP{j},FRE{j},AMP_all,'linear');
    DAM_int(:,j)=interp1(AMP{j},DAM{j},AMP_all,'linear');
end

%% calculate the averaged and 95% CI for frequency and damping ratio
AMP_avg=AMP_all;
for loop1=1:size(FRE_int,1)
    FRE_avg(loop1,1)=mean(FRE_int(loop1,:));
    FRE_std(loop1,1)=std(FRE_int(loop1,:));
    FRE_upper(loop1,1)=FRE_avg(loop1,1)+1.96*FRE_std(loop1,1);
    FRE_lower(loop1,1)=FRE_avg(loop1,1)-1.96*FRE_std(loop1,1);
    
    DAM_avg(loop1,1)=mean(DAM_int(loop1,:));
    DAM_std(loop1,1)=std(DAM_int(loop1,:));
    DAM_upper(loop1,1)=DAM_avg(loop1,1)+1.96*DAM_std(loop1,1);
    DAM_lower(loop1,1)=DAM_avg(loop1,1)-1.96*DAM_std(loop1,1);
end

end

