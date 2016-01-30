% function [h,chi2Stat,pvalue]=chi2Fun(mean1,var1,mean2,var2)
% 
% Computes a chi squared test given the means and SD, assuming unequal variances
%
% INPUT:  mean1, mean2 = array of means (of length N)
%         var1, var2   = array of variances (e.g. variance=diag(Cov))
%
% OUTPUT: h      = test result: 1 if different, 0 if same, NaN if invalid 
%         chi2   = array of chi2 statistics
%         pvalue = array of p-values (if p<0.05/N the result is
%         significant)
% 
% Luca Mazzucato January 2016

function [h,chi2Stat,pvalue]=chi2Fun(mean1,var1,mean2,var2)

if size(mean1,2)>size(mean1,1)
    mean1=mean1';
end
if size(var1,2)>size(var1,1)
    var1=var1';
end
if size(mean2,2)>size(mean2,1)
    mean2=mean2';
end
if size(var2,2)>size(var2,1)
    var2=var2';
end

% invalid entries
meanElim=isnan(mean1) | isnan(mean2);
varianceElim=(isnan(var1) | var1<=0) | (isnan(var2) | var2<=0);
indElim=meanElim | varianceElim;

NumElim=sum(indElim);
if NumElim>0
    fprintf('ttestFun: %d entries have zero variance or are NaNs...\n',NumElim);    
    fprintf('          corresponding outputs set to NaN\n');
end

chi2Stat=(mean1-mean2).^2./(var1+var2);
pvalue=1-chi2cdf(chi2Stat,2);
h=pvalue<0.05/numel(mean1);
pvalue(indElim)=NaN;
