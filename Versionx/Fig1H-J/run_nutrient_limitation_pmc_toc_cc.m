function [growthRate_PMC,ppGpp_PMC,...
          growthRate_TOC,ppGpp_TOC,...
          growthRate_CC, ppGpp_CC] = ...
          run_nutrient_limitation_pmc_toc_cc(hp,nutrList,ppGppFixed)

%   hp: host cell parameters
%   nutrList: list of nutrient values

%   set Matlab solvers
tol = 1e-6;
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

growthRate_PMC  = zeros(length(nutrList),1);
ppGpp_PMC       = zeros(length(nutrList),1);
initCond        = zeros(length(nutrList),3);    %   used as initial guess for TOC and CC model

%   initial condition
x0 = [10,10,10]; 
tic;
[~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(end),0,hp);
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
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC,x0,options_fsolve,nutrList(i),0,hp);
        if (exitflag <=0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(i),0,hp);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr(i),0,hp);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    x0 = x;
    ppGpp_PMC(i)  = x0(3);
    initCond(i,:) = x0;
    
    [~,growthRate_PMC(i)] = Ecoli_GR_ODE_PMC(0,x0,nutrList(i),0,hp);
    if (growthRate_PMC(i)<tol)
        break;
    end
end

%%  theoretical optimal control
growthRate_TOC          = zeros(length(nutrList),1);
ppGpp_TOC               = zeros(length(nutrList),1);

%   we run the simulation from high nutrient to low nutrient and stop it
%   when the growth rate becomes zero
for k=length(nutrList):-1:1
    k
    [optPpGpp,~,exitflag] = fmincon(@optFunc,initCond(k,3),[],[],[],[],[0],[1000],[],options_fmincon,nutrList(k),0,hp,initCond(k,[1,2]));
    assert(exitflag>0);

    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,initCond(k,[1,2]),options_ode15s_CC,nutrList(k),0,hp,optPpGpp);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    ppGpp_TOC(k) = optPpGpp;
    [~,growthRate_TOC(k)] = Ecoli_GR_ODE_CC(0,x(end,:),nutrList(k),0,hp,optPpGpp);
    if (growthRate_TOC(k)<tol)
        break;
    end
end

%%  constant control

growthRate_CC = zeros(length(nutrList),1);
ppGpp_CC      = zeros(length(nutrList),1);

%   we run the simulation from high nutrient to low nutrient and stop it
%   when the growth rate becomes zero
for i=length(nutrList):-1:1
    i
    %   choose to use ode15s or fsolve
    if (1)
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC,initCond(i,[1,2]),options_fsolve,nutrList(i),0,hp,ppGppFixed);
        if (exitflag <=0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,initCond(i,[1,2]),options_ode15s_CC,nutrList(i),0,hp,ppGppFixed);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,initCond(i,[1,2]),options_ode15s_CC,nutrList(i),0,hp,ppGppFixed);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    
    ppGpp_CC(i) = ppGppFixed;
    [~,growthRate_CC(i)] = Ecoli_GR_ODE_CC(0,x,nutrList(i),0,hp,ppGppFixed);
    if (growthRate_CC(i)<tol)
        ppGpp_CC(1:i) = ppGppFixed;
        break;
    end
end

end

%   maximization of growth rate
function growthRate = optFunc(ppGppFixed,nutr,cm,hp,initCond)

tol   = 1e-6;
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
