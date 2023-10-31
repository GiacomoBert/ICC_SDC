%% Compute the Smallest Detectable Change (SDC) as in
% Schambra HM,OgdenRT, Martínez-Hernández IE,LinX,Chang YB, Rahman A, Edwards DJ and Krakauer JW(2015)
%The reliability of repeated TMS measures in older adults and inpatients with subacute  and chronicstroke.
% Front. Cell.Neurosci.9:335. doi: 10.3389/fncel.2015.00335
%
% SDC_out=SDC(x1,x2,n)
% inputs:
% x1 and x2 must be column vectors
% n is a scalar indicating the number of subjects
% outputs:
% SDC_out=SDC
% nSDC_out=normalized SDC by dividing it to the grandmean
% SEMas=standard error of the measurament


% SDC
% function [SEMas, SDC, nSDC]=SDC(x1,x2,n)
function [SDC_out, nSDC_out, SEMas_out]=SDC(x1,x2,n)

%% SDC amp
data_T1=x1;
data_T2=x2;
subnum=n;

data_T1_T2=[data_T1 data_T2]; %full matrix
grand_mean=mean(data_T1_T2, 'all'); %grand mean
mean_data_T1=mean(data_T1); %conditional mean T1
mean_data_T2=mean(data_T2); %conditional mean T2
% total error sum|(X-grandmean)^2
for count=1:numel(data_T1_T2)
    dummy_total(count)=(data_T1_T2(count)-grand_mean)^2;
end
SS_total=sum(dummy_total, 'all');

% between subj error k*sum|(subj_mean-grandeman)^2
df_between=subnum-1;
subj_means=(mean(data_T1_T2, 2));
for count=1:length(subj_means)
    dummy_between(count)=(subj_means(count)-grand_mean)^2;
end
SS_between=sum(dummy_between, 'all')*2;

% within subj error (trials(days) error - residual error)
%trials error (sum| subjnum *(mean_condition - grandmean)^2)
df_trials=2-1; %num conditions-1
dummy_trials_1=subnum*((mean_data_T1-grand_mean)^2);
dummy_trials_2=subnum*((mean_data_T2-grand_mean)^2);
SS_trials=dummy_trials_1+dummy_trials_2;

%residual error (total error - between error - trials error)
SS_residual=SS_total-SS_between-SS_trials;
df_residual=df_between-df_trials;

% within subj error
SS_within=SS_residual+SS_trials;
df_within=df_trials+df_residual;

%mean sum square
MS_between=SS_between/df_between;
MS_within=SS_within/df_within;
MS_trials=SS_trials/df_trials;
MS_residual=SS_residual/df_residual;

%SEMas standard error of the measurament
SEMas_out=sqrt(MS_within);
SDC_out=SEMas_out*sqrt(2)*1.96; %SDC corrected for two way design and at 95%
nSDC_out=(SDC_out/abs(grand_mean)); % normalization;



