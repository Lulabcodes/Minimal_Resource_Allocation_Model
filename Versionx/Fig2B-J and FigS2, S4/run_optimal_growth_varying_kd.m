function [growthRate,growthRate_opt] = run_optimal_growth_varying_kd(hp,nutrStar,kd_fc,eta,gamma)

%   hp: host cell parameters

growthRate          = zeros(4,length(kd_fc));
growthRate_opt      = zeros(4,length(kd_fc));

%   set Matlab solvers
tol = 1e-10;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC_unsatVer);
options_fmincon= optimset('Display','none','TolFun',tol,'TolX',tol,'Algorithm','active-set');

%   time span
tspan = [0 10^10];

for Flag=0:3
 
    %   initial condition
    x0 = [10,10,10];
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_unsatVer,tspan,x0,options_ode15s_PMC,nutrStar,0,hp,Flag,eta,gamma);
    if (isempty(te))
        x0 = x(end,:);
    else
        error('Error: Oscillation Detected for the current parameter set!');
    end
    
    hp_copy = hp;
    for i=1:length(kd_fc)
        [Flag,i]
        
        hp_copy.('kdPpGpp') = hp.('kdPpGpp')*kd_fc(i);
        
        %%  simulation using full model
        if (0)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC_unsatVer,x0,options_fsolve,nutrStar,0,hp_copy,Flag,eta,gamma);
            if (exitflag <=0 || sum(x<=0) > 0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_unsatVer,tspan,x0,options_ode15s_PMC,nutrStar,0,hp_copy,Flag,eta,gamma);
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_unsatVer,tspan,x0,options_ode15s_PMC,nutrStar,0,hp_copy,Flag,eta,gamma);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        x0 = x;
        
        [~,growthRate(Flag+1,i)] = Ecoli_GR_ODE_PMC_unsatVer(0,x0,nutrStar,0,hp_copy,Flag,eta,gamma);
        
        %%  optimal control model
        [optPpGpp,~,exitflag] = fmincon(@optFunc,x0(3),[],[],[],[],[0],[1500],[],options_fmincon,nutrStar,0,hp_copy,x0([1,2]),Flag,eta,gamma);
        if (exitflag<0)
            options_fmincon_alter= optimset('Display','none','TolFun',tol2,'TolX',tol2,'Algorithm','interior-point');
            [optPpGpp,~,exitflag] = fmincon(@optFunc,x0(3),[],[],[],[],[0],[1500],[],options_fmincon_alter,nutrStar,0,hp_copy,x0([1,2]),Flag,eta,gamma);
            assert(exitflag>0);
        end
        
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,x0([1,2]),options_ode15s_CC,nutrStar,0,hp_copy,optPpGpp,Flag,eta,gamma);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        [~,growthRate_opt(Flag+1,i)] = Ecoli_GR_ODE_CC_unsatVer(0,x(end,:),nutrStar,0,hp_copy,optPpGpp,Flag,eta,gamma);
    end
end

end

%   maximization of growth rate
function growthRate = optFunc(ppGppFixed,nutr,cm,hp,initCond,Flag,eta,gamma)

tol   = 1e-10;
tspan = [0 10^10];
option_ode15s_CC = odeset('NonNegative',[1,2],'RelTol',tol,'AbsTol',tol,'Events',@myEvent_CC_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

if (1)
    [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_unsatVer,initCond,options_fsolve,nutr,cm,hp,ppGppFixed,Flag,eta,gamma);
    if (exitflag <=0 || sum(x<=0) > 0)
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,initCond,option_ode15s_CC,nutr,cm,hp,ppGppFixed,Flag,eta,gamma);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
else
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,initCond,option_ode15s_CC,nutr,cm,hp,ppGppFixed,Flag,eta,gamma);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    x = x(end,:);   %   only keep the steady state solution
end
[~,growthRate] = Ecoli_GR_ODE_CC_unsatVer(0,x,nutr,cm,hp,ppGppFixed,Flag,eta,gamma);

growthRate = -growthRate;

end

