function [AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate] = ...
    run_cell_growth_vary_ppGpp_diff_ka(hp,nutrStar,ppGppList,ppGppStar,initCond,kaToShow)

%   hp: host cell parameters
%   nutrStar: the specific nutrient level
%   ppGppList: values of fixed ppGpp levels

growthRate          = zeros(length(kaToShow),length(ppGppList));
AminoAcid           = zeros(length(kaToShow),length(ppGppList));
Ribosome            = zeros(length(kaToShow),length(ppGppList));
ppGppSynRate        = zeros(length(kaToShow),length(ppGppList));
ppGppDegRate        = zeros(length(kaToShow),length(ppGppList));

%   set Matlab solvers
tol = 1e-10;
options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol,'TolFun',tol);

%   time span
tspan = [0 10^10];

%   nearest index of ppGppStar
[~,IndexPpGppStar] = min(abs(ppGppList-ppGppStar));

%   running rightward
hp_copy = hp;
for k=1:length(kaToShow)
    
    k
    
    hp_copy.('thetaAaAc') = kaToShow(k);
    
    %   initial condition
    if (k==1)
        x0 = initCond;
    else
        x0 = [AminoAcid(k-1,IndexPpGppStar),Ribosome(k-1,IndexPpGppStar)]
    end
    
    for i=IndexPpGppStar:length(ppGppList)
        
        %   choose to use ode15s or fsolve
        if (1)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC,x0,options_fsolve,nutrStar,0,hp_copy,ppGppList(i));
            if (exitflag <=0 || sum(x<=0) > 0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,nutrStar,0,hp_copy,ppGppList(i));
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,nutrStar,0,hp_copy,ppGppList(i));
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        AminoAcid(k,i)        = x(1);
        Ribosome(k,i)         = x(2);
        [~,growthRate(k,i),~,~,~,ppGppSynRate(k,i),ppGppDegRate(k,i)] = ...
            Ecoli_GR_ODE_CC(0,x,nutrStar,0,hp_copy,ppGppList(i));
        if (growthRate(k,i) < tol)
            growthRate(k,i) = 0;
        end
        x0 = x;
    end
end

%   running leftward
hp_copy = hp;
for k=1:length(kaToShow)
    
    k
    
    hp_copy.('thetaAaAc') = kaToShow(k);
    
    %   initial condition
    if (k==1)
        x0 = initCond;
    else
        x0 = [AminoAcid(k-1,IndexPpGppStar),Ribosome(k-1,IndexPpGppStar)]
    end
    
    for i=IndexPpGppStar:-1:1

        %   choose to use ode15s or fsolve
        if (1)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_CC,x0,options_fsolve,nutrStar,0,hp_copy,ppGppList(i));
            if (exitflag <=0 || sum(x<=0) > 0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,nutrStar,0,hp_copy,ppGppList(i));
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,nutrStar,0,hp_copy,ppGppList(i));
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        AminoAcid(k,i)        = x(1);
        Ribosome(k,i)         = x(2);
        [~,growthRate(k,i),~,~,~,ppGppSynRate(k,i),ppGppDegRate(k,i)] = ...
            Ecoli_GR_ODE_CC(0,x,nutrStar,0,hp_copy,ppGppList(i));
        if (growthRate(k,i) < tol)
            growthRate(k,i) = 0;
        end
        x0 = x;
    end
end

end

