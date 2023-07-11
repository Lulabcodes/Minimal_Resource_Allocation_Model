clear;
clc;

addpath('../');

%%  read parameters
hostparams = readParameters();

%%  Use different models for ppGpp-dependent ribosome synthesis function

nutrList = [0,10.^[0:0.1:3]];
[growthRate_PMC_model,ppGpp_PMC_model,growthRate_TOC_model,ppGpp_TOC_model] = ...
    run_nutrient_limitation_pmc_toc_diff_model(hostparams,nutrList);

%%  plot growth rate optimization for different models of RCF
%   Fig2K

figure();
hold on;

%   maker size
MS = 8;

%   line width
LW = 1.5;

%   different model
model = {'hill2';'linear';'exponential';'parabolic'};

%   color code
ColorPalette = [255,109,109;255,219,107;72,255,167;61,191,255;198,73,255]/255;

for q=1:length(model)
    
    indexPMC_nonzero = find(growthRate_PMC_model(q,:),1,'first');
    indexTOC_nonzero = find(growthRate_TOC_model(q,:),1,'first');
    
    subplot(2,2,q);
    hold on;
    
    plot(nutrList,growthRate_PMC_model(q,:),'k-','LineWidth',LW','Color',ColorPalette(q,:));
    plot(nutrList,growthRate_TOC_model(q,:),'k--','LineWidth',LW','Color',ColorPalette(q,:));
    legend({'PMC','TOC'},'Location','southeast');
    
    axis square;
    box on;
    xlabel('Nutrient uptake rate (h^{-1})');
    ylabel('Growth rate (h^{-1})');
    axis([0,1000,0,2.4]);
    set(gca,'YTick',[0,1.2,2.4]);
    set(gca,'XTick',[0,500,1000]);
    
    title(model{q});
end
