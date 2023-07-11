function [res_PMC,res_TOC] = run_growth_rate_optimization_vary_Ka(hp,nutrStar,kaList)

%   hp: host cell parameters

%   set Matlab solvers
tol = 1e-10;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
                            'RelTol',tol,...
                            'AbsTol',tol,...
                            'Events',@myEvent_PMC);
options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC);
options_fmincon= optimset('Display','off','TolFun',tol,'TolX',tol);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span                        
tspan = [0 10^10];

%%  ppGpp-mediated control
%   res_PMC are arranged in the order of
%   AminoAcid_PMC, Ribosome_PMC, ppGpp_PMC, growthRate_PMC

res_PMC  = zeros(length(kaList),6);

%   initial condition
x0 = [10,10,10]; 
hp_copy = hp;
hp_copy.('thetaAaAc') = kaList(1);
tic;
[~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrStar,0,hp_copy);
if (isempty(te))
    x0 = x(end,:);
else
    error('Error: Oscillation Detected for the current parameter set!');
end


for i=1:length(kaList)   
    i

    hp_copy.('thetaAaAc') = kaList(i);
    
    %   choose to use ode15s or fsolve
    if (1)
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC,x0,options_fsolve,nutrStar,0,hp_copy);
        if (exitflag <=0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrStar,0,hp_copy);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrStar,0,hp_copy);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    x0 = x;
    res_PMC(i,1:3)= x0;
    
    [~,res_PMC(i,4)] = Ecoli_GR_ODE_PMC(0,x0,nutrStar,0,hp_copy);
end

%%  theoretical optimal control
res_TOC  = zeros(length(kaList),6);

xguess_ppGpp = res_PMC(1,3);
xguess0      = res_PMC(1,[1,2]);

for k=1:length(kaList)
    k
    
    hp_copy.('thetaAaAc') = kaList(k);

    [optPpGpp,~,exitflag] = fmincon(@optFunc,xguess_ppGpp,[],[],[],[],[0],[1000],[],options_fmincon,nutrStar,0,hp_copy,xguess0);
    assert(exitflag>0);
    xguess_ppGpp = optPpGpp;
    res_TOC(k,3) = optPpGpp;

    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,xguess0,options_ode15s_CC,nutrStar,0,hp_copy,optPpGpp);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    xguess0 = x(end,:);
    res_TOC(k,[1,2]) = xguess0;
    
    [~,res_TOC(k,4)] = Ecoli_GR_ODE_CC(0,x(end,:),nutrStar,0,hp_copy,optPpGpp);
end

end

%   maximization of growth rate
function growthRate = optFunc(ppGppFixed,nutr,cm,hp,initCond)

tol   = 1e-10;
tspan = [0 10^10];
option_ode15s_CC = odeset('NonNegative',[1,2],'RelTol',tol,'AbsTol',tol,'Events',@myEvent_CC);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

[x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC,initCond,options_fsolve,nutr,cm,hp,ppGppFixed);
if (exitflag <=0 || sum(x<=0)>0)
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,initCond,option_ode15s_CC,nutr,cm,hp,ppGppFixed);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    x = x(end,:);   %   only keep the steady state solution
end
[~,growthRate] = Ecoli_GR_ALE_CC(x,nutr,cm,hp,ppGppFixed);

growthRate = -growthRate;

end
