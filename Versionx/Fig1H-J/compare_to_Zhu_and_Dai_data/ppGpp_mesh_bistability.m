clear;
clc;

addpath('../');
addpath('../../');

%%  read parameters
hostparams = readParameters();

%   conversion factor between OD600 and um^3
convFac     = 0.50;

% Zhu and Dai 2019 NAR
ZD2019_Growth_relA      = [1.725, 1.210, 0.856, 0.240]; % 1/h
ZD2019_ppGpp_relA       = [16.602, 26.380, 49.434, 147.620]/convFac;
ZD2019_ppGpp_relA_err   = [0, 0, 15.202, 17.054]/convFac; 

ZD2019_Growth_mesh      = [0.609, 0.407, 0.281, 0.179]; % 1/h
ZD2019_ppGpp_mesh       = [74.539,45.241,26.695,14.872]/convFac;
ZD2019_ppGpp_mesh_err   = [12.5380,9.5570,9.2760,6.9450]/convFac; 

figure();
hold on;

desired_growthRate = [0.609, 0.609];
parameters = {'k_s';'d_g'};

for q=1:2
    
    %%  find nutrient levels
    [nutr,growthRate] = run_nutrient_limitation_pmc(hostparams);
    [~,I] = min(abs(growthRate-desired_growthRate(q)));
    nutrStar = nutr(I);

    if (q==1)
        paraRatio  = 10.^[0:0.005:5];
    else
        paraRatio  = 10.^[5:-0.005:0];
    end
    [AminoAcid_PMC,Ribosome_PMC,ppGpp_PMC,growthRate_PMC,stability_PMC] = ...
            run_one_parameter_scan_PMC_ks_kd(hostparams,nutrStar,paraRatio);

    hold on;
    plot(ppGpp_PMC(2,:),growthRate_PMC(2,:),'k-');

    axis square;
    box on;
    xlabel('ppGpp concentration (\muM)');
    ylabel('Growth rate (1/h)');
    %axis([0.1,10,0,110]);
    %set(gca,'XTick',[0.1,1,10]);
    %set(gca,'YTick',[0,50,100]);
    title(parameters{q});
    set(gca,'XMinorTick','off');
    set(gca,'YMinorTick','off');
end

errorbar(ZD2019_ppGpp_mesh, ZD2019_Growth_mesh, [], [], ZD2019_ppGpp_mesh_err, ZD2019_ppGpp_mesh_err, 'ko');
