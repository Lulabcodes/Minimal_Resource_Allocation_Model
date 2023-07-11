clear;
clc;

addpath('../');

%%  read parameters
hostparams = readParameters();

%% find nutrient level that correspondes to growth rate 1.0
[nutr,growthRate,AminoAcid,Ribosome,ppGpp] = run_nutrient_limitation(hostparams);
[growthRate,I] = unique(growthRate);
nutr           = nutr(I);
AminoAcid      = AminoAcid(I);
Ribosome       = Ribosome(I);
ppGpp          = ppGpp(I);

nutrStar        = pchip(growthRate,nutr,1.0);
AminoAcidStar   = pchip(nutr,AminoAcid,nutrStar);
RibosomeStar    = pchip(nutr,Ribosome,nutrStar);
ppGppStar       = pchip(nutr,ppGpp,nutrStar);
initCond        = [AminoAcidStar,RibosomeStar];

%%  run cell growth by varying ppGpp

ppGppList = [0,10.^[1:0.001:3]];
eta     = 0.1;
gamma   = 0.1;
[AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth] = ...
    run_cell_growth_vary_ppGpp(hostparams,nutrStar,ppGppList,ppGppStar,initCond,eta,gamma);

%   save file
save('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

%%  plot growth rate, ribosome, amino acids and ppGpp synth./deg. rate at varied ppGpp concentrations
%   Fig2B,C (the first column)
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');
plot_cell_growth_vary_ppGpp(ppGppList,AminoAcid,Ribosome,growthRate,...
    ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth,hostparams);

%%  plot ultrasensitivity values at varied ppGpp concentrations
%   Fig2D
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');
plot_sensitivity_values_vary_ppGpp(ppGppList,AminoAcid,Ribosome,growthRate,...
    ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth,hostparams,0);

%%  run rate-balance analysis using full model
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

ppGppToShow       = [40:10:110];
AminoAcidToShow   = 10.^[-2:0.01:5];

[AminoAcidSynRate,AminoAcidDegRate] = ...
    run_rate_balance_analysis_full_model(hostparams,nutrStar,ppGppToShow,AminoAcidToShow,eta,gamma);

%   save to file
save('data_rate_balance_analysis_full_model_gr1.0_eta0.1_gamma0.1.mat');

%%  plot rate-balance analysis using full model
%   Fig 2F-H
load('data_rate_balance_analysis_full_model_gr1.0_eta0.1_gamma0.1.mat');

%   ppGppList and AminoAcid are computed from previous steps
plot_rate_balance_analysis_full_model(...
    ppGppToShow,AminoAcidToShow,AminoAcidSynRate,AminoAcidDegRate,ppGppList,AminoAcid,IndexMaxGrowth);

%%  show robustness of growth rate optimization at high ultrasensitivity
%   Fig2E

load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');
plot_optimal_growth_varying_kd(hostparams,ppGppList,growthRate,ppGppSynRate,IndexMaxGrowth,IndexLowerLimit,IndexUpperLimit);

%%  vary ppGpp-mediated feedback types
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

ppGppList_fb = [10.^[-3:0.01:5]];
[AminoAcid_fb,Ribosome_fb,growthRate_fb,ppGppSynRate_fb,ppGppDegRate_fb,...
    IndexLowerLimit_fb,IndexUpperLimit_fb,IndexMaxGrowth_fb] = ...
    run_ppGpp_rate_balance_diff_feedback(hostparams,nutrStar,ppGppList_fb,ppGppStar,initCond);

save('data_growth_rate_different_fb_type_gr1.0.mat');

%%  plot ppGpp reaction rates and growth rates using different feedback types
%   Fig 2J
load('data_growth_rate_different_fb_type_gr1.0.mat');

plot_rate_balance_different_feedback_type(...
    ppGppList_fb,ppGppSynRate_fb,ppGppDegRate_fb,growthRate_fb,IndexLowerLimit_fb,IndexUpperLimit_fb);

%%  vary ka to show the effects of ultrasensitivity on growth rate optimization
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

kaList = 10.^[-1:0.1:5,5.01:0.01:6,6.001:0.001:6.5,6.50001:0.0001:7];
[res_PMC,res_TOC] = run_growth_rate_optimization_vary_Ka(hostparams,nutrStar,kaList);

save('data_growth_rate_optimization_vary_Ka_gr1.0.mat');

% %%  plot growth rate optimization at varied Ka values
% load('data_growth_rate_optimization_vary_Ka_gr1.0.mat');
% plot_growth_rate_optimization_vary_Ka(kaList,res_PMC,res_TOC);

%%  run to find relationship between Ka and ultrasensitivity (it works until 10^5.6
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

kaScan = 10.^[0:0.1:3];
[AminoAcid_kascan,Ribosome_kascan,growthRate_kascan,ppGppSynRate_kascan,ppGppDegRate_kascan] = ...
    run_cell_growth_vary_ppGpp_diff_ka(hostparams,nutrStar,ppGppList,ppGppStar,initCond,kaScan);

save('data_ultrasensitivity_at_diff_ka_gr1.0_10power5.6.mat');

%%  plot relationship between Ka and ultrasensitivity
%   Fig2I
load('data_growth_rate_optimization_vary_Ka_gr1.0.mat');
load('data_ultrasensitivity_at_diff_ka_gr1.0.mat')

plot_function_ka_ultrasensitivity(kaScan,ppGppList,ppGppSynRate_kascan,growthRate_kascan,...
    kaList,res_PMC,res_TOC);

%%  run ppGpp variation at large amino acid toxicity

ppGppList = [0,10.^[1:0.001:3]];
eta     = 0.1;
gamma   = 0.1;
hostparams.thetaAaTox = 1e+2;

[AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth] = ...
    run_cell_growth_vary_ppGpp(hostparams,nutrStar,ppGppList,ppGppStar,initCond,eta,gamma);

%   save file
save('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1_AAtox1e2.mat');

%%  run rate-balance analysis using full model
%   amino acid toxicity is large
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1_AAtox1e2.mat');
hostparams.thetaAaTox = 1e+2;

ppGppToShow       = [40:10:110];
AminoAcidToShow   = 10.^[-2:0.01:5];

[AminoAcidSynRate,AminoAcidDegRate] = ...
    run_rate_balance_analysis_full_model(hostparams,nutrStar,ppGppToShow,AminoAcidToShow,eta,gamma);

%   save to file
save('data_rate_balance_analysis_full_model_gr1.0_eta0.1_gamma0.1_AAtox1e2.mat');

%%  plot rate-balance analysis using full model
%   FigS2A
load('data_rate_balance_analysis_full_model_gr1.0_eta0.1_gamma0.1_AAtox1e2.mat');

%   ppGppList and AminoAcid are computed from previous steps
plot_rate_balance_analysis_full_model_vary_Ka(...
    ppGppToShow,AminoAcidToShow,AminoAcidSynRate,AminoAcidDegRate,ppGppList,AminoAcid,IndexMaxGrowth);

%   overlay with wild type curve
load('data_rate_balance_analysis_full_model_gr1.0_eta0.1_gamma0.1.mat');
Flag = 0;
deriv               = fnder(spline(ppGppList,AminoAcid(Flag+1,:)),1);
localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./AminoAcid(Flag+1,:));
plot(ppGppList,localSensitivity,'k-');
plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,40],'k--');

%%  show comparison of cell growth at different Ka
load('data_cell_growth_vary_ppGpp_gr1.0_eta0.1_gamma0.1.mat');

kaToShow = 10.^[1,2,3,4];
[AminoAcid_ka,Ribosome_ka,growthRate_ka,ppGppSynRate_ka,ppGppDegRate_ka,...
    IndexLowerLimit_ka,IndexUpperLimit_ka,IndexMaxGrowth_ka] = ...
    run_cell_growth_vary_ppGpp_diff_ka(hostparams,nutrStar,ppGppList,ppGppStar,initCond,kaToShow);

save('data_cell_growth_vary_ppGpp_diff_ka_gr1.0.mat');

%%  plot comparison of cell growth at different Ka
%   FigS2B
load('data_cell_growth_vary_ppGpp_diff_ka_gr1.0.mat');

plot_cell_growth_vary_ppGpp_diff_Ka(kaToShow,ppGppList,...
    AminoAcid_ka,Ribosome_ka,growthRate_ka,ppGppSynRate_ka,ppGppDegRate_ka,...
    IndexLowerLimit_ka,IndexUpperLimit_ka,IndexMaxGrowth_ka,hostparams);

%%  calculate band width and ultrasensitivity at different nutrient quality

selected_GR    = [0.1:0.05:1.9];
BW             = zeros(length(selected_GR),1);
US             = zeros(length(selected_GR),1);
ppGppList      = [0,10.^[1:0.001:3]];

for i=1:length(selected_GR)
    i
    
    [nutr,growthRate,AminoAcid,Ribosome,ppGpp] = run_nutrient_limitation(hostparams);
    [growthRate,I] = unique(growthRate);
    nutr           = nutr(I);
    AminoAcid      = AminoAcid(I);
    Ribosome       = Ribosome(I);
    ppGpp          = ppGpp(I);
    
    nutrStar        = pchip(growthRate,nutr,selected_GR(i));
    AminoAcidStar   = pchip(nutr,AminoAcid,nutrStar);
    RibosomeStar    = pchip(nutr,Ribosome,nutrStar);
    ppGppStar       = pchip(nutr,ppGpp,nutrStar);
    initCond        = [AminoAcidStar,RibosomeStar];
    
    [AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth] = ...
        run_cell_growth_vary_ppGpp_bw(hostparams,nutrStar,ppGppList,ppGppStar,initCond,0,0);
    
    %   calculate bandwidth
    [maxGR,I] = max(growthRate);
    [~,indexL] = min(abs(growthRate(1:I)-maxGR/2));
    [~,indexR] = min(abs(growthRate(I+1:length(growthRate))-maxGR/2));
    indexR = indexR+I;
    BW(i) = log10(AminoAcid(indexR)/AminoAcid(indexL));
    
    %   calculate ultrasensitivity
    deriv    = fnder(spline(ppGppList,AminoAcid),1);
    localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./AminoAcid);
    US(i) = localSensitivity(I);
    
end

save('data_relat_bandw_ultras_vary_nutr.mat');

%%  plot relationship between bandwidth and ultrasensitivity
%   FigS4B,C
load('data_relat_bandw_ultras_vary_nutr.mat');

figure();

subplot(1,2,2);
hold on;

plot(selected_GR,BW./log10(US),'ko-','MarkerFaceColor','r');
axis square;
box on;
xlabel('Growth rate (h^{-1})');
ylabel('B_w/log_{10}(S_a^*)');
axis([0,2,2,3]);
set(gca,'XTick',[0,1,2]);
set(gca,'YTick',[2,2.5,3]);

subplot(1,2,1);
hold on;

[nutr,growthRate,AminoAcid,Ribosome,ppGpp] = run_nutrient_limitation(hostparams);
[growthRate,I] = unique(growthRate);
nutr           = nutr(I);
AminoAcid      = AminoAcid(I);
Ribosome       = Ribosome(I);
ppGpp          = ppGpp(I);
[~,index1p5_ss]   = min(abs(growthRate-1.5));
aa_ss = AminoAcid(index1p5_ss);
gr_ss = growthRate(index1p5_ss);
plot(aa_ss,gr_ss,'ko','MarkerFaceColor','r');

nutrStar        = pchip(growthRate,nutr,1.5);
AminoAcidStar   = pchip(nutr,AminoAcid,nutrStar);
RibosomeStar    = pchip(nutr,Ribosome,nutrStar);
ppGppStar       = pchip(nutr,ppGpp,nutrStar);
initCond        = [AminoAcidStar,RibosomeStar];

[AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth] = ...
    run_cell_growth_vary_ppGpp_bw(hostparams,nutrStar,ppGppList,ppGppStar,initCond,0,0);

AA_smooth = 10.^[-2:0.01:6];
GR_smooth = pchip(AminoAcid,growthRate,AA_smooth);
plot(AA_smooth,GR_smooth,'k-');
hold on;
plot([1e-2,1e6],[max(growthRate),max(growthRate)],'k--');
axis square;
box on;
axis([1e-2,1e+6,0,1.6]);
set(gca,'XScale','log');
set(gca,'XTick',[1e-2,1e2,1e6]);
set(gca,'XTicklabel',{'10^{-2}';'10^{2}';'10^6'});
set(gca,'YTick',[0,0.8,1.6]);
xlabel('Amino acid (\muM)');
ylabel('Growth rate (h^{-1})');

