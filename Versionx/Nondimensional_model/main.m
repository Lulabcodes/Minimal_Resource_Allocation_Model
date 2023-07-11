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

desiredGR      = 1.0;
nutrStar        = pchip(growthRate,nutr,desiredGR);
AminoAcidStar   = pchip(nutr,AminoAcid,nutrStar);
RibosomeStar    = pchip(nutr,Ribosome,nutrStar);
ppGppStar       = pchip(nutr,ppGpp,nutrStar);
initCond        = [AminoAcidStar,RibosomeStar];

%%  check result
normR = RibosomeStar*hostparams.massR/hostparams.beta;
normA = AminoAcidStar/hostparams.beta;
normGR = hostparams.massR*desiredGR/hostparams.kelongMax;
normke = nutrStar*hostparams.massR/hostparams.kelongMax/hostparams.massE;
normKa = hostparams.thetaAaAc/hostparams.beta;
normKc = hostparams.KdAcTRNA*hostparams.massR/hostparams.ratioTR/hostparams.beta;
normKu = hostparams.KdUAcTRNA*hostparams.massR/hostparams.ratioTR/hostparams.beta;
normks = hostparams.ksPpGpp/hostparams.kelongMax;

normR_pred = 2*normke*normKc*(normA+normKa)*0.5/(...
    (normke*normKc*(normA+normKa)-normke*normA*0.5-normke*normKa*normKc*0.5/normKu)+...
    sqrt((normke*normKc*(normA+normKa)-normke*normA*0.5-normke*normKa*normKc*0.5/normKu)^2+...
    4*normke*normKc*0.5*(normA+normKa)*(0.05*normA+normA*normA)))


    
%%  run cell growth by varying ppGpp

ppGppList = [0,10.^[1:0.01:3]];
[AminoAcid,Ribosome,growthRate,ppGppSynRate,ppGppDegRate] = ...
    run_cell_growth_vary_ppGpp(hostparams,nutrStar,ppGppList,ppGppStar,initCond,0,0);

%%  plot

normR = Ribosome(1,:)*hostparams.massR/hostparams.beta;
normA = AminoAcid(1,:)/hostparams.beta;
normGR = hostparams.massR*growthRate(1,:)/hostparams.kelongMax;
normke = nutrStar*hostparams.massR/hostparams.kelongMax/hostparams.massE;
normKa = hostparams.thetaAaAc/hostparams.beta;
normKc = hostparams.KdAcTRNA*hostparams.massR/hostparams.ratioTR/hostparams.beta;
normKu = hostparams.KdUAcTRNA*hostparams.massR/hostparams.ratioTR/hostparams.beta;
normks = hostparams.ksPpGpp/hostparams.kelongMax;
ar     = hostparams.fAa+normke*(1+normKc/hostparams.phiRMax);
al     = normke*normKa*normKc*(1/normKu+1/hostparams.phiRMax)./ar;
lambdaub = normke/ar;

%   ppGpp sensitivity
subplot(1,3,1);
deriv               = fnder(spline(ppGppList,Ribosome(1,:)),1);
localS              = abs(ppval(deriv,ppGppList).*ppGppList./Ribosome(1,:));
localS_theory       = normR/hostparams.phiRMax;

plot(ppGppList,localS,'k-');
hold on;
plot(ppGppList,localS_theory,'k--');
xlim([50,150]);

subplot(1,3,2);

plot(AminoAcid(1,:),growthRate(1,:));
set(gca,'XScale','log');
xlim([1e-5,1e+5]);

%   amino acid sensitivity
subplot(1,3,3);
deriv               = fnder(spline(ppGppList,AminoAcid(1,:)),1);
localS              = abs(ppval(deriv,ppGppList).*ppGppList./AminoAcid(1,:));
localS_theory = normR./(1-normR)./sqrt((ar.*(1-normR)-1).^2-4*al.*ar.*(1-normR).^2).*normR/hostparams.phiRMax;

plot(ppGppList,localS,'k-');
hold on;
plot(ppGppList,localS_theory,'k--');
xlim([50,150]);
