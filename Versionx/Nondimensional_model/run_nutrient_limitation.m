function [nutr,growthRate,AminoAcid,Ribosome,ppGpp,kelong,fracCharge,fracRactive] = ...
          run_nutrient_limitation(hp)

%   hp: host cell parameters

nutr                = 10.^[-4:0.05:4];
AminoAcid           = zeros(length(nutr),1);
Ribosome            = zeros(length(nutr),1);
ppGpp               = zeros(length(nutr),1);
growthRate          = zeros(length(nutr),1);
kelong              = zeros(length(nutr),1);
fracCharge          = zeros(length(nutr),1);
fracRactive         = zeros(length(nutr),1);

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
                            'RelTol',tol,...
                            'AbsTol',tol,...
                            'Events',@myEvent_PMC);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span                        
tspan = [0 10^10];

%   initial condition
x0 = [10,10,10]; 
tic;
[~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr(end),0,hp);
if (isempty(te))
    x0 = x(end,:);
else
    error('Error: Oscillation Detected for the current parameter set!');
end

for i=length(nutr):-1:1
    %   choose to use ode15s or fsolve
    if (1)
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC,x0,options_fsolve,nutr(i),0,hp);
        if (exitflag <=0 || sum(x<=0)>0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr(i),0,hp);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr(i),0,hp);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    x0 = x;
    
    AminoAcid(i)    = x0(1);
    Ribosome(i)     = x0(2);
    ppGpp(i)        = x0(3);
    [~,growthRate(i),kelong(i),fracCharge(i),fracRactive(i)] = Ecoli_GR_ODE_PMC(0,x0,nutr(i),0,hp);
    if (growthRate(i)<tol)
        break;
    end
end

end

