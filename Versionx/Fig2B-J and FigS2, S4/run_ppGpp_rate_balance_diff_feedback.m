function [AminoAcid_fb,Ribosome_fb,growthRate_fb,ppGppSynRate_fb,ppGppDegRate_fb,...
    IndexLowerLimit_fb,IndexUpperLimit_fb,IndexMaxGrowth_fb] = ...
    run_ppGpp_rate_balance_diff_feedback(hp,nutrStar,ppGppList,ppGppStar,initCond)

%   hp: host cell parameters
%   nutrStar: the specific nutrient level
%   ppGppList: values of fixed ppGpp levels

growthRate_fb          = zeros(3,length(ppGppList));
AminoAcid_fb           = zeros(3,length(ppGppList));
Ribosome_fb            = zeros(3,length(ppGppList));
ppGppSynRate_fb        = zeros(3,length(ppGppList));
ppGppDegRate_fb        = zeros(3,length(ppGppList));

%   set Matlab solvers
tol = 1e-10;
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
for Flag = 0:2
    
    %   initial condition
    x0 = initCond;
    
    for i=IndexPpGppStar:length(ppGppList)
        
        [Flag,i]
        %   choose to use ode15s or fsolve
        if (1)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_Feedback,x0,options_fsolve,nutrStar,0,hp,ppGppList(i),Flag);
            if (exitflag <=0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Feedback,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag);
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Feedback,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        AminoAcid_fb(Flag+1,i)        = x(1);
        Ribosome_fb(Flag+1,i)         = x(2);
        [~,growthRate_fb(Flag+1,i),~,~,~,ppGppSynRate_fb(Flag+1,i),ppGppDegRate_fb(Flag+1,i)] = ...
            Ecoli_GR_ODE_CC_Feedback(0,x,nutrStar,0,hp,ppGppList(i),Flag);
        if (growthRate_fb(Flag+1,i) < tol)
            growthRate_fb(Flag+1,i) = 0;
        end
        x0 = x;
    end
end

%   running leftward
for Flag=0:2
    %   initial condition
    x0 = initCond;
    
    for i=IndexPpGppStar:-1:1
        
        [Flag,i]
        
        %   choose to use ode15s or fsolve
        if (0)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC_Feedback,x0,options_fsolve,nutrStar,0,hp,ppGppList(i),Flag);
            if (exitflag <=0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Feedback,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag);
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC_Feedback,tspan,x0,options_ode15s_CC,nutrStar,0,hp,ppGppList(i),Flag);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        AminoAcid_fb(Flag+1,i)        = x(1);
        Ribosome_fb(Flag+1,i)         = x(2);
        [~,growthRate_fb(Flag+1,i),~,~,~,ppGppSynRate_fb(Flag+1,i),ppGppDegRate_fb(Flag+1,i)] = ...
            Ecoli_GR_ODE_CC_Feedback(0,x,nutrStar,0,hp,ppGppList(i),Flag);
        if (growthRate_fb(Flag+1,i) < tol)
            growthRate_fb(Flag+1,i) = 0;
        end
        x0 = x;
    end
end

%  Find index of ppGpp for maximum growth rate, lower and upper limits
IndexMaxGrowth_fb  = zeros(1,3);
IndexLowerLimit_fb = zeros(1,3);
IndexUpperLimit_fb = zeros(1,3);
for Flag=0:2
    [~,IndexMaxGrowth_fb(Flag+1)]  = max(growthRate_fb(Flag+1,:));
    IndexLowerLimit_fb(Flag+1)     = find(growthRate_fb(Flag+1,1:IndexMaxGrowth_fb(Flag+1)),1,'first');
    IndexUpperLimit_fb(Flag+1)     = find(growthRate_fb(Flag+1,IndexMaxGrowth_fb(Flag+1):end),1,'last')+IndexMaxGrowth_fb(Flag+1)-1;
end

end

