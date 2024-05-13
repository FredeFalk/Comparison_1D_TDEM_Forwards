function out = get_analytic_tdem_circLoop(models,S,LoopArea)
    
% S. H. Ward and G. W. Hohmann. Electromagnetic theory for geophysical
% applications. In M. N. Nabighian, editor, Electromagnetic Methods in
% Applied Geophysics, Volume 1, Theory, pages 131–311. Society of
% Exploration Geophysicists, 1987.

% Concentric loop transmitter of radius r on homogenous halfspace with
% conductivity sigma_h. dbdt z component.

% Waveform is a step-off waveform

Nmod = numel(models);
r = sqrt(LoopArea/pi);
%gates = 

mu0 = 1.2566E-6;

constant_factor = @(sigma,r) 1./(1*sigma*power(r,3));

tau = @(t,sigma,r) r*sqrt((sigma*mu0)./(4*t));

time_dependent_factor = @(t,sigma,r) 3*erf(tau(t,sigma,r))-(2/sqrt(pi)).*tau(t,sigma,r).*(3+2.*power(tau(t,sigma,r),2)).*exp(-power(tau(t,sigma,r),2));

GateArray = S.General.GateArray;
gates = GateArray(:,1);

NGate = numel(gates);

dBdt = zeros(Nmod,NGate);

for i = 1:Nmod
    curmod = models{i};
    cond = 1./curmod(2,1);

    cur_dbdt = constant_factor(cond,r)*time_dependent_factor(gates,cond,r);
    dBdt(i,:) = cur_dbdt/(LoopArea);
end

%OUTPUT: dBdt array has dimensions G-by-N where G is the number of gates
dBdt_LM = dBdt;
dBdt_HM = dBdt;

out.HM.dBdt = dBdt_HM;
out.HM.Gates = gates;
out.LM.dBdt = dBdt_LM;
out.LM.Gates = gates;
out.LM.UseGates = gates*0+1;
out.HM.UseGates = gates*0+1;

end