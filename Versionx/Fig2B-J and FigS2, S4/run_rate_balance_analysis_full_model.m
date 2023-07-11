function [AminoAcidSynRate,AminoAcidDegRate] = ...
    run_rate_balance_analysis_full_model(hp,nutrStar,ppGppToShow,AminoAcidToShow,eta,gamma)

%   set Matlab solvers
tol = 1e-8;
options_ode15s_PMC = odeset('NonNegative',[1],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_GA_unsatVer);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span
tspan = [0 10^10];

AminoAcidSynRate = zeros(4,length(ppGppToShow),length(AminoAcidToShow));
AminoAcidDegRate = zeros(4,length(ppGppToShow),length(AminoAcidToShow));

for Flag=0:3
    
    for k=1:length(ppGppToShow)
        
        [Flag,k]
        
        %   initial condition
        found = 1;
        while (found)
            x0 = rand*hp.beta*hp.phiRMax/hp.massR;
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_GA_unsatVer,tspan,x0,options_ode15s_PMC,...
                nutrStar,0,hp,Flag,eta,gamma,ppGppToShow(k),AminoAcidToShow(end));
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x0 = x(end,:);
            if (x0>0)
                found = 0;
            end
        end

        for j=length(AminoAcidToShow):-1:1
            if (0)
                [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_GA_unsatVer,x0,options_fsolve,...
                    nutrStar,0,hp,Flag,eta,gamma,ppGppToShow(k),AminoAcidToShow(j));
                if (exitflag <=0 || sum(x<=0)>0)
                    tic;
                    [~,x,te] = ode15s(@Ecoli_GR_ODE_GA_unsatVer,tspan,x0,options_ode15s_PMC,...
                        nutrStar,0,hp,Flag,eta,gamma,ppGppToShow(k),AminoAcidToShow(j));
                    if (~isempty(te))
                        error('Error: Oscillation Detected for the current parameter set!');
                    end
                    x = x(end,:);   %   only keep the steady state solution
                end
            else
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_GA_unsatVer,tspan,x0,options_ode15s_PMC,...
                    nutrStar,0,hp,Flag,eta,gamma,ppGppToShow(k),AminoAcidToShow(j));
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
            x0 = x;
            [~,~,AminoAcidSynRate(Flag+1,k,j),AminoAcidDegRate(Flag+1,k,j)] = ...
                Ecoli_GR_ODE_GA_unsatVer(0,x0,nutrStar,0,hp,Flag,eta,gamma,ppGppToShow(k),AminoAcidToShow(j));
        end
        
    end
    
end

end

