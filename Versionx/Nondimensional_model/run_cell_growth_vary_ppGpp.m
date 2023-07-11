function [AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate] = ...
    run_cell_growth_vary_ppGpp(hp,nutrStar,ppGppList,ppGppStar,initCond,eta,gamma)

%   hp: host cell parameters
%   nutrStar: the specific nutrient level
%   ppGppList: values of fixed ppGpp levels

growthRate          = zeros(1,length(ppGppList));
AminoAcid           = zeros(1,length(ppGppList));
Ribosome            = zeros(1,length(ppGppList));
ppGppSynRate        = zeros(1,length(ppGppList));
ppGppDegRate        = zeros(1,length(ppGppList));

%   set Matlab solvers
tol = 1e-8;
options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span
tspan = [0 10^10];

Flag = 0;

%   nearest index of ppGppStar
[~,IndexPpGppStar] = min(abs(ppGppList-ppGppStar));

%   initial condition
x0 = initCond;

for i=IndexPpGppStar:length(ppGppList)
    i
    %   choose to use ode15s or fsolve
    if (0)
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_unsatVer,x0,options_fsolve,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
        if (exitflag <=0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    AminoAcid(Flag+1,i)        = x(1);
    Ribosome(Flag+1,i)         = x(2);
    [~,growthRate(Flag+1,i),~,~,~,ppGppSynRate(Flag+1,i),ppGppDegRate(Flag+1,i)] = ...
        Ecoli_GR_ODE_CC_unsatVer(0,x,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
    if (growthRate(Flag+1,i) < tol)
        growthRate(Flag+1,i) = 0;
    end
    x0 = x;
end


%   initial condition
x0 = initCond;

for i=IndexPpGppStar:-1:1
    
    [Flag,i]
    
    %   choose to use ode15s or fsolve
    if (0)
        [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_unsatVer,x0,options_fsolve,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
        if (exitflag <=0)
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
    else
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_unsatVer,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
    end
    AminoAcid(Flag+1,i)        = x(1);
    Ribosome(Flag+1,i)         = x(2);
    [~,growthRate(Flag+1,i),~,~,~,ppGppSynRate(Flag+1,i),ppGppDegRate(Flag+1,i)] = ...
        Ecoli_GR_ODE_CC_unsatVer(0,x,nutrStar,0,hp,ppGppList(i),Flag,eta,gamma);
    if (growthRate(Flag+1,i) < tol)
        growthRate(Flag+1,i) = 0;
    end
    x0 = x;
end



end

