clear;
clc;

addpath('../');

%%  read parameters
hostparams = readParameters();

%%  run nutrient limitation for pmc, toc, and cc at steady state

%   Important!!!
%   The initial condition for upshift/downshit is selected where the growth
%   rates of pmc, toc and cc are equal
%   This particular point can be found through very dense sampling of nutrient quality
%   For example, nutrList = 10.^[-1:0.01:1.89,1.90:0.0001:2,2.1:0.01:3];
%   After doing this, we found that the particular spot is nutr = 92.02,
%   corresponding to 10^(1.9639), which must be included in the sampling
nutrList = 10.^[-1:0.1:1.9,1.92:0.01:1.96,1.9639,1.97:0.01:1.99,2:0.1:3];

ppGppFixed = 60;   %   uM
[growthRate_PMC_ss,ppGpp_PMC_ss,...
    growthRate_TOC_ss,ppGpp_TOC_ss,...
    growthRate_CC_ss, ppGpp_CC_ss] = ...
    run_nutrient_limitation_pmc_toc_cc(hostparams,nutrList,ppGppFixed);

%%  run nutrient upshift/downshift

[growthRate_PMC_ss_unique,I] = unique(growthRate_PMC_ss);
nutrList_unique = nutrList(I);
growthRate_TOC_ss_unique = growthRate_TOC_ss(I);
growthRate_CC_ss_unique  = growthRate_CC_ss(I);

%   upshift:    shiftNutr(1) -> critNutr
%   downshift:  critNutr -> shiftNutr(2)
[~,J]           = min(abs(growthRate_TOC_ss_unique(2:end)-growthRate_CC_ss_unique(2:end)));
critNutr        = nutrList_unique(J);
shiftNutr       = pchip(growthRate_PMC_ss_unique,nutrList_unique,[0.8,1.4]);
shiftNutr       = [shiftNutr,critNutr];
shiftPpGpp_TOC  = pchip(nutrList,ppGpp_TOC_ss,shiftNutr);
shiftPpGpp_CC   = pchip(nutrList,ppGpp_CC_ss,shiftNutr);

[t_PMC_upshift,growthRate_PMC_t_upshift,ppGpp_PMC_t_upshift,...
    t_PMC_downshift,growthRate_PMC_t_downshift,ppGpp_PMC_t_downshift,...
    t_TOC_upshift,growthRate_TOC_t_upshift,ppGpp_TOC_t_upshift,...
    t_TOC_downshift,growthRate_TOC_t_downshift,ppGpp_TOC_t_downshift,...
    t_CC_upshift,growthRate_CC_t_upshift,ppGpp_CC_t_upshift,...
    t_CC_downshift,growthRate_CC_t_downshift,ppGpp_CC_t_downshift] ...
    = run_nutrient_limitation_upshift_downshift(hostparams,shiftNutr,shiftPpGpp_TOC,shiftPpGpp_CC);

%%   plot Fig1H-J

figure();

%   color code
ColorPalette = [72,255,167;61,191,255;198,73,255]/255;

%   maker size
MS = 8;

%   line width
LW = 1.5;

%   steady state growth rate
subplot(2,3,1);
hold on;

indexPMC_nonzero = find(growthRate_PMC_ss,1,'first');
indexTOC_nonzero = find(growthRate_TOC_ss,1,'first');
indexCC_nonzero  = find(growthRate_CC_ss,1,'first');

plot(nutrList,growthRate_PMC_ss,'k-','LineWidth',LW','Color',ColorPalette(1,:));
plot(nutrList,growthRate_TOC_ss,'k-','LineWidth',LW','Color',ColorPalette(2,:));
plot(nutrList(1:indexCC_nonzero-1),growthRate_CC_ss(1:indexCC_nonzero-1),'k-','LineWidth',LW','Color',ColorPalette(3,:));
plot(nutrList(indexCC_nonzero:end),growthRate_CC_ss(indexCC_nonzero:end),'k-','LineWidth',LW','Color',ColorPalette(3,:));
xcorr_middle = (nutrList(indexCC_nonzero-1)+nutrList(indexCC_nonzero))/2;
plot([xcorr_middle,xcorr_middle],[0,growthRate_CC_ss(indexCC_nonzero)],'k--','LineWidth',LW','Color',ColorPalette(3,:));
%legend({'PMC';'TOC';'CC'},'Location','southeast');

%   nutrient upshift
plot([shiftNutr(1),shiftNutr(1)],[0,2.4],'k--');
plot([shiftNutr(2),shiftNutr(2)],[0,2.4],'k--');
plot([critNutr(1),critNutr(1)],[0,2.4],'k--');

axis square;
box on;
xlabel('Nutrient quality');
ylabel('Growth rate (h^{-1})');
axis([0,300,0,2]);
set(gca,'YTick',[0,1,2]);
set(gca,'XTick',[0:100:300]);

%   steady state growth rate
subplot(2,3,4);
hold on;

plot(nutrList(indexPMC_nonzero:end),hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_PMC_ss(indexPMC_nonzero:end)),'k-','LineWidth',LW','Color',ColorPalette(1,:));
plot(nutrList(indexTOC_nonzero:end),hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_TOC_ss(indexTOC_nonzero:end)),'k-','LineWidth',LW','Color',ColorPalette(2,:));
plot(nutrList, hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_CC_ss),'k-','LineWidth',LW','Color',ColorPalette(3,:));
%legend({'PMC';'TOC';'CC'},'Location','southeast');

axis square;
box on;
xlabel('Nutrient quality');
ylabel('RCF');
axis([0,300,0,1]);
set(gca,'YTick',[0,0.5,1]);
set(gca,'XTick',[0:100:300]);

%   growth rate dynamics: upshift
subplot(2,3,2);
hold on;

plot(t_PMC_upshift,growthRate_PMC_t_upshift,'k-','LineWidth',LW,'Color',ColorPalette(1,:));
plot(t_TOC_upshift,growthRate_TOC_t_upshift,'k-','LineWidth',LW,'Color',ColorPalette(2,:));
plot(t_CC_upshift, growthRate_CC_t_upshift, 'k-','LineWidth',LW,'Color',ColorPalette(3,:));
plot([-0.4,4],[1.199,1.199],'k--');

axis square;
box on;
axis([-0.4,4,1.1,1.5]);
xlabel('Time (h)');
ylabel('Growth rate (h^{-1})');
set(gca,'XTick',[0,2,4]);
set(gca,'YTick',[1.1,1.3,1.5]);

%   ppGpp dynamics: upshift
subplot(2,3,5);
hold on;

plot(t_PMC_upshift,hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_PMC_t_upshift),'k-','LineWidth',LW,'Color',ColorPalette(1,:));
plot(t_TOC_upshift,hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_TOC_t_upshift),'k-','LineWidth',LW,'Color',ColorPalette(2,:));
plot(t_CC_upshift, hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_CC_t_upshift), 'k-','LineWidth',LW,'Color',ColorPalette(3,:));

axis square;
box on;
axis([-0.4,4,0.4,1.0]);
xlabel('Time (h)');
ylabel('RCF');
set(gca,'XTick',[0,2,4]);
set(gca,'YTick',[0.4,0.7,1.0]);

%   growth rate dynamics: downshift
subplot(2,3,3);
hold on;

plot(t_PMC_downshift,growthRate_PMC_t_downshift,'k-','LineWidth',LW,'Color',ColorPalette(1,:));
plot(t_TOC_downshift,growthRate_TOC_t_downshift,'k-','LineWidth',LW,'Color',ColorPalette(2,:));
plot(t_CC_downshift, growthRate_CC_t_downshift, 'k-','LineWidth',LW,'Color',ColorPalette(3,:));

axis square;
box on;
axis([-0.4,4,0.4,1.2]);
xlabel('Time (h)');
ylabel('Growth rate (h^{-1})');
set(gca,'XTick',[0,2,4]);
set(gca,'YTick',[0.4,0.8,1.2]);

%   ppGpp dynamics: downshift
subplot(2,3,6);
hold on;

plot(t_PMC_downshift,hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_PMC_t_downshift),'k-','LineWidth',LW,'Color',ColorPalette(1,:));
plot(t_TOC_downshift,hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_TOC_t_downshift),'k-','LineWidth',LW,'Color',ColorPalette(2,:));
plot(t_CC_downshift, hostparams.('thetaPpGppR')./(hostparams.('thetaPpGppR')+ppGpp_CC_t_downshift), 'k-','LineWidth',LW,'Color',ColorPalette(3,:));

axis square;
box on;
axis([-0.4,4,0,0.6]);
xlabel('Time (h)');
ylabel('RCF');
set(gca,'XTick',[0,2,4]);
set(gca,'YTick',[0,0.3,0.6]);

% %%  Use different models for ppGpp-dependent ribosome synthesis function
% 
% %nutrList                = 10.^[-1:0.01:1.89,1.90:0.0001:2,2.1:0.01:3];
% nutrList = [0,10.^[0:0.1:3]];
% [growthRate_PMC_model,ppGpp_PMC_model,growthRate_TOC_model,ppGpp_TOC_model] = ...
%     run_nutrient_limitation_pmc_toc_diff_model(hostparams,nutrList);
% 
% % save to file
% save('steady_state_pmc_toc_nutrient_limitation_diff_model.mat');
% 
% %%  plot growth rate optimization for different models of RCF
% 
% load('steady_state_pmc_toc_nutrient_limitation_diff_model.mat');
% 
% figure();
% hold on;
% 
% %   maker size
% MS = 8;
% 
% %   line width
% LW = 1.5;
% 
% %   different model
% model = {'hill2';'linear';'exponential';'parabolic'};
% 
% %   color code
% ColorPalette = [255,109,109;255,219,107;72,255,167;61,191,255;198,73,255]/255;
% 
% for q=1:length(model)
%     
%     indexPMC_nonzero = find(growthRate_PMC_model(q,:),1,'first');
%     indexTOC_nonzero = find(growthRate_TOC_model(q,:),1,'first');
%     
%     subplot(2,2,q);
%     hold on;
%     
%     plot(nutrList,growthRate_PMC_model(q,:),'k-','LineWidth',LW','Color',ColorPalette(q,:));
%     plot(nutrList,growthRate_TOC_model(q,:),'k--','LineWidth',LW','Color',ColorPalette(q,:));
%     legend({'PMC','TOC'},'Location','southeast');
%     
%     axis square;
%     box on;
%     xlabel('Nutrient uptake rate (h^{-1})');
%     ylabel('Growth rate (h^{-1})');
%     axis([0,1000,0,2.4]);
%     set(gca,'YTick',[0,1.2,2.4]);
%     set(gca,'XTick',[0,500,1000]);
%     
%     title(model{q});
% end
