%% Compute the Smallest Detectable Change (SDC) as in
% Schambra HM,OgdenRT, Martínez-Hernández IE,LinX,Chang YB, Rahman A, Edwards DJ and Krakauer JW(2015)
%The reliability of repeated TMS measures in older adults and inpatients with subacute  and chronicstroke.
% Front. Cell.Neurosci.9:335. doi: 10.3389/fncel.2015.00335
%
% SDC_out=SDC(x1,x2,n)
% inputs:
% x1 and x2 must be column vectors and contains the data of the two session that need to be compared
% each row is a subject. 
% n is a scalar indicating the number of subjects
% nboot is the numebr of iteration for the bootstrapped CI
% outputs:
% SDC_out=SDC
% nSDC_out=normalized SDC by dividing it to the grandmean
% SEMas=standard error of the measurament
% SDC_CI_95= upper and lower bound of the SDC confidence interval at 95% bootstapped

function [SDC_out, nSDC_out, SEMas_out, SDC_CI_95, nSDC_CI_95]=SDC_CI(x1,x2,n,nboot)

[SDC_out, nSDC_out, SEMas_out]=SDC(x1,x2,n); %compute SDC 

%compute SDC CI 
mySDC = @(x1,x2,n) SDC(x1,x2,n);
stat_set=statset('UseParallel',1);
SDC_CI_95=bootci(nboot,{mySDC,x1,x2,n},'type', 'norm', 'Options', stat_set);
nSDC_CI_95=SDC_CI_95/abs(mean([x1 x2], 'all'));


    