function [growthRate_PMC,ppGpp_PMC,growthRate_TOC,ppGpp_TOC] = ...
    run_nutrient_limitation_pmc_toc_diff_model(hp,nutrList)

%   hp: host cell parameters
%   nutrList: list of nutrient values

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC_Rcf = odeset('NonNegative',[1,2,3],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC_Rcf);
options_ode15s_CC_Rcf = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC_Rcf);
options_fmincon= optimset('Display','off','TolFun',tol,'TolX',tol);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span
tspan = [0 10^10];

%   model
model = {'hill2';'linear';'exponential';'parabolic';'wt'};

%%  ppGpp-mediated control

growthRate_PMC  = zeros(length(model),length(nutrList));
ppGpp_PMC       = zeros(length(model),length(nutrList));
initCond        = zeros(length(model),length(nutrList),3);    %   used as initial guess for TOC and CC model

for q=1:length(model)
    
    %   initial condition
    x0 = [10,10,10];
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_Rcf,tspan,x0,options_ode15s_PMC_Rcf,nutrList(end),0,hp,model{q});
    if (isempty(te))
        x0 = x(end,:);
    else
        error('Error: Oscillation Detected for the current parameter set!');
    end
    
    %   we run the simulation from high nutrient to low nutrient and stop it
    %   when the growth rate becomes zero
    for i=length(nutrList):-1:1
        i
        %   choose to use ode15s or fsolve
        if (1)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC_Rcf,x0,options_fsolve,nutrList(i),0,hp,model{q});
            if (exitflag <=0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_Rcf,tspan,x0,options_ode15s_PMC_Rcf,nutrList(i),0,hp,model{q});
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC_Rcf,tspan,x0,options_ode15s_PMC,nutr(i),0,hp,model{q});
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        x0 = x;
        ppGpp_PMC(q,i)  = x0(3);
        initCond(q,i,:) = x0;
        
        [~,growthRate_PMC(q,i)] = Ecoli_GR_ODE_PMC_Rcf(0,x0,nutrList(i),0,hp,model{q});
        if (growthRate_PMC(q,i)<tol)
            break;
        end
    end
    
end

%%  theoretical optimal control
growthRate_TOC          = zeros(length(model),length(nutrList));
ppGpp_TOC               = zeros(length(model),length(nutrList));

for q=1:length(model)
    
    %   we run the simulation from high nutrient to low nutrient and stop it
    %   when the growth rate becomes zero
    for k=length(nutrList):-1:1
        k
        [optPpGpp,~,exitflag] = fmincon(@optFunc,initCond(q,k,3),[],[],[],[],[0],[1000],[],options_fmincon,nutrList(k),0,hp,initCond(q,k,[1,2]),model{q});
        assert(exitflag>0);
        
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Rcf,tspan,initCond(q,k,[1,2]),options_ode15s_CC_Rcf,nutrList(k),0,hp,optPpGpp,model{q});
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        ppGpp_TOC(q,k) = optPpGpp;
        [~,growthRate_TOC(q,k)] = Ecoli_GR_ODE_CC_Rcf(0,x(end,:),nutrList(k),0,hp,optPpGpp,model{q});
        if (growthRate_TOC(q,k)<tol)
            break;
        end
    end
end

end

%   maximization of growth rate
function growthRate = optFunc(ppGppFixed,nutr,cm,hp,initCond,model)

tol   = 1e-6;
tspan = [0 10^10];
option_ode15s_CC = odeset('NonNegative',[1,2],'RelTol',tol,'AbsTol',tol,'Events',@myEvent_CC_Rcf);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

[x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_Rcf,initCond,options_fsolve,nutr,cm,hp,ppGppFixed,model);
if (exitflag <=0 || sum(x<=0)>0)
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Rcf,tspan,initCond,option_ode15s_CC,nutr,cm,hp,ppGppFixed,model);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    x = x(end,:);   %   only keep the steady state solution
end
[~,growthRate] = Ecoli_GR_ALE_CC_Rcf(x,nutr,cm,hp,ppGppFixed,model);

growthRate = -growthRate;

end
