function [relerrors,residuals,maxsigma] = get_relerror(Reference,M,UseGatesRef,UseGatesM)
%Sigma is 3% and maxsigma is assuming that SigmaA=3% and SigmaB=3%
sigma = 3;

%Usegates is a vector that should have the same length both for the
%reference and M. The code starts by determining which gates to compare
%between Reference and M.
Nmods = size(Reference,1);

for i = 1:Nmods
   cur_ref = Reference(i,:);
   cur_mod = M(i,:);
   
   VecA = double(UseGatesRef*0);
   VecB = VecA;

   VecA(UseGatesRef>0) = cur_ref;
   VecB(UseGatesM>0) = cur_mod;

   cur_residual = VecB-VecA;
   cur_relerror = cur_residual./VecA;

   residuals(:,i) = cur_residual;
   relerrors(:,i) = cur_relerror;
end
%%
MeansB = VecA*0;
MeansA = VecA*0;

MeansA(UseGatesM>0) = mean(M);
MeansB(UseGatesRef>0) = mean(Reference);

fraction = power(MeansA/MeansB,2);
maxsigma = sqrt(sigma^2*(1+fraction));
%%
end