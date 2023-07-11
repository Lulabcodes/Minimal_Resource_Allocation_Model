function [AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth] = ...
    run_cell_growth_vary_ppGpp(hp,nutrStar,ppGppList,ppGppStar,initCond,eta,gamma)

%   hp: host cell parameters
%   nutrStar: the specific nutrient level
%   ppGppList: values of fixed ppGpp levels

growthRate          = zeros(4,length(ppGppList));
AminoAcid           = zeros(4,length(ppGppList));
Ribosome            = zeros(4,length(ppGppList));
ppGppSynRate        = zeros(4,length(ppGppList));
ppGppDegRate        = zeros(4,length(ppGppList));

%   set Matlab solvers
tol = 1e-8;
options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span
tspan = [0 10^10];

%   nearest index of ppGppStar
[~,IndexPpGppStar] = min(abs(ppGppList-ppGppStar));

%   running rightward
for Flag = 0:3
    
    %   initial condition
    x0 = initCond;
    
    for i=IndexPpGppStar:length(ppGppList)
        
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

%   running leftward
for Flag=0:3
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

%  Find index of ppGpp for maximum growth rate, lower and upper limits
IndexMaxGrowth  = zeros(1,4);
IndexLowerLimit = zeros(1,4);
IndexUpperLimit = zeros(1,4);
for Flag=0:3
    [~,IndexMaxGrowth(Flag+1)]  = max(growthRate(Flag+1,:));
    IndexLowerLimit(Flag+1)     = find(growthRate(Flag+1,1:IndexMaxGrowth(Flag+1)),1,'first');
    IndexUpperLimit(Flag+1)     = find(growthRate(Flag+1,IndexMaxGrowth(Flag+1):end),1,'last')+IndexMaxGrowth(Flag+1)-1;
end

end

