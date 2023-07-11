function [cm,growthRate,AminoAcid,Ribosome,ppGpp,kelong,fracCharge,fracRactive] = ...
    run_translation_limitation(hp,nutrList)

assert(length(nutrList)==4);

%   nutrList: list of nutrient values
%   hp: host cell parameters

cm                  = [0:0.2:5,6:1:48,50:0.5:100];
AminoAcid           = zeros(length(nutrList),length(cm));
Ribosome            = zeros(length(nutrList),length(cm));
ppGpp               = zeros(length(nutrList),length(cm));
growthRate          = zeros(length(nutrList),length(cm));
kelong              = zeros(length(nutrList),length(cm));
fracCharge          = zeros(length(nutrList),length(cm));
fracRactive         = zeros(length(nutrList),length(cm));

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC);
options_fsolve = optimoptions('fsolve','Display','none','TolX',tol);

%   time span
tspan = [0 10^10];

indexCut = ones(length(cm),1)*length(cm);
for k=1:length(nutrList)
    
    %   initial condition
    x0 = [10,10,10];
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(k),0,hp);
    if (isempty(te))
        x0 = x(end,:);
    else
        error('Error: Oscillation Detected for the current parameter set!');
    end
    
    %   we run the simulation from high nutrient to low nutrient and stop it
    %   when the growth rate becomes zero
    for i=1:length(cm)
        
        %   choose to use ode15s or fsolve
        if (1)
            [x,~,exitflag] = fsolve(@Ecoli_GR_ALE_PMC,x0,options_fsolve,nutrList(k),cm(i),hp);
            if (exitflag <=0 || sum(x<=0)>0)
                tic;
                [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(k),cm(i),hp);
                if (~isempty(te))
                    error('Error: Oscillation Detected for the current parameter set!');
                end
                x = x(end,:);   %   only keep the steady state solution
            end
        else
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(k),cm(i),hp);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
        end
        x0 = x;
        
        AminoAcid(k,i)    = x0(1);
        Ribosome(k,i)     = x0(2);
        ppGpp(k,i)        = x0(3);
        [~,growthRate(k,i),kelong(k,i),fracCharge(k,i),fracRactive(k,i)] = Ecoli_GR_ODE_PMC(0,x0,nutrList(k),cm(i),hp);
        if (growthRate(k,i)<tol)
            indexCut(k) = i-1;
            break;
        end
    end
end

%   R-protein fraction (affiliated proteins excluded)
phiR = Ribosome*hp.massR./hp.beta/1.6;

%   plot
figure();

%   All data are extracted from the following paper:
%   Paper: Reduction of translating ribosomes enables Escherichi coli to
%   maintain elongation rates during slow growth
%   Xiongfeng Dai et al. (2017)

%   color code
ColorPalette = [255,109,109;255,219,107;72,255,167;61,191,255;198,73,255]/255;

%   maker size
MS = 8;

%   line width
LW = 1.5;

%----------------------------------------------------------------
%  R protein fraction (affliated protein excluded) vs Growth rate
%----------------------------------------------------------------

subplot(1,2,1);
hold on;

for k=1:length(nutrList)
    plot(growthRate(k,1:indexCut(k)),phiR(k,1:indexCut(k)),'k-','LineWidth',LW);
end
axis square;
box on;
xlabel('Growth rate (hr^{-1})');
ylabel('R fraction');

%   RDM+0.2%Glucose
growthRate_RDM = [1.8,1.46,1.08,0.87,0.57];
phiR_RDM       = [0.476,0.495,0.551,0.582,0.621]*0.76/1.6;
plot(growthRate_RDM,phiR_RDM,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(1,:));
RDM_cm = [0,2,4,6,8];
growthRate_Pred = pchip(cm(1:indexCut(1)),growthRate(1,1:indexCut(1)),RDM_cm);
phiR_Pred = pchip(cm(1:indexCut(1)),phiR(1,1:indexCut(1)),RDM_cm);
plot(growthRate_Pred,phiR_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(1,:));

%   0.2% Glucose + 10mM NH4Cl
growthRate_Glu  = [0.98,0.71,0.53,0.41,0.33,0.26];
phiR_Glu        = [0.294,0.358,0.440,0.487,0.511,0.569]*0.76/1.6;
plot(growthRate_Glu,phiR_Glu,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(2,:));
Glu_cm = [0,2,4,6,8,9];
growthRate_Pred = pchip(cm(1:indexCut(2)),growthRate(2,1:indexCut(2)),Glu_cm);
phiR_Pred = pchip(cm(1:indexCut(2)),phiR(2,1:indexCut(2)),Glu_cm);
plot(growthRate_Pred,phiR_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(2,:));

%   0.2% Fructose + 10 mM NH4Cl
growthRate_Fru  = [0.69,0.46,0.35,0.27,0.21];
phiR_Fru        = [0.217,0.286,0.323,0.394,0.457]*0.76/1.6;
plot(growthRate_Fru,phiR_Fru,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(3,:));
Fru_cm = [0,2,4,6,8];
growthRate_Pred = pchip(cm(1:indexCut(3)),growthRate(3,1:indexCut(3)),Fru_cm);
phiR_Pred = pchip(cm(1:indexCut(3)),phiR(3,1:indexCut(3)),Fru_cm);
plot(growthRate_Pred,phiR_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(3,:));

%   60 mM acetate + 10mM NH4Cl
growthRate_Ace  = [0.46,0.25,0.18];
phiR_Ace        = [0.172,0.246,0.304]*0.76/1.6;
plot(growthRate_Ace,phiR_Ace,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(4,:));
Ace_cm = [0,3,6];
growthRate_Pred = pchip(cm(1:indexCut(4)),growthRate(4,1:indexCut(4)),Ace_cm);
phiR_Pred = pchip(cm(1:indexCut(4)),phiR(4,1:indexCut(4)),Ace_cm);
plot(growthRate_Pred,phiR_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(4,:));

xlim([0,2]);
ylim([0,0.4]);

%----------------------------------------------------------------
%   Peptide elongation rate vs Growth rate
%----------------------------------------------------------------

subplot(1,2,2);
hold on;
for k=1:length(nutrList)
    plot(growthRate(k,1:indexCut(k)),kelong(k,1:indexCut(k))/3600,'k-','LineWidth',LW);
end
axis square;
box on
xlabel('Growth rate (hr^{-1})');
ylabel('Peptide synthesis rate (aa/s)');

%   RDM+0.2%Glucose
growthRate_RDM = [1.8,1.08,0.57];
kelong_RDM     = [16.7,16.8,17.3];
plot(growthRate_RDM,kelong_RDM,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(1,:));
RDM_cm = [0,4,8];
growthRate_Pred = pchip(cm(1:indexCut(1)),growthRate(1,1:indexCut(1)),RDM_cm);
kelong_Pred = pchip(cm(1:indexCut(1)),kelong(1,1:indexCut(1))/3600,RDM_cm);
plot(growthRate_Pred,kelong_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(1,:));

%   0.2% Glucose + 10mM NH4Cl
growthRate_Glu = [0.98,0.71,0.53,0.41,0.33,0.26];
kelong_Glu     = [15.9,16.0,16.1,16.2,16.5,16.6];
plot(growthRate_Glu,kelong_Glu,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(2,:));
Glu_cm = [0,2,4,6,8,9];
growthRate_Pred = pchip(cm(1:indexCut(2)),growthRate(2,1:indexCut(2)),Glu_cm);
kelong_Pred = pchip(cm(1:indexCut(2)),kelong(2,1:indexCut(2))/3600,Glu_cm);
plot(growthRate_Pred,kelong_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(2,:));

%   0.2% Fructose + 10 mM NH4Cl
growthRate_Fru = [0.69,0.35,0.21];
kelong_Glu     = [14.7,15.9,16.3];
plot(growthRate_Fru,kelong_Glu,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(3,:));
Fru_cm = [0,4,8];
growthRate_Pred = pchip(cm(1:indexCut(3)),growthRate(3,1:indexCut(3)),Fru_cm);
kelong_Pred = pchip(cm(1:indexCut(3)),kelong(3,1:indexCut(3))/3600,Fru_cm);
plot(growthRate_Pred,kelong_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(3,:));

%   60 mM acetate + 10mM NH4Cl
growthRate_Ace = [0.46,0.25,0.17];
kelong_Ace     = [12.6,14.5,15.6];
plot(growthRate_Ace,kelong_Ace,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(4,:));
Ace_cm = [0,3,6];
growthRate_Pred = pchip(cm(1:indexCut(4)),growthRate(4,1:indexCut(4)),Ace_cm);
kelong_Pred = pchip(cm(1:indexCut(4)),kelong(4,1:indexCut(4))/3600,Ace_cm);
plot(growthRate_Pred,kelong_Pred,'kx','MarkerSize',MS,'LineWidth',LW,'MarkerEdgeColor',ColorPalette(4,:));

xlim([0,2]);
ylim([10,20]);

end

