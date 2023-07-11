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

desired_growthRate = [1.725, 0.609];
parameters = {'k_s';'d_g'};

for q=1:2
    
    %%  find nutrient levels
    [nutr,growthRate] = run_nutrient_limitation_pmc(hostparams);
    [~,I] = min(abs(growthRate-desired_growthRate(q)));
    nutrStar = nutr(I);

    paraRatio  = 10.^[0:0.005:5];
    [AminoAcid_PMC,Ribosome_PMC,ppGpp_PMC,growthRate_PMC,stability_PMC] = ...
            run_one_parameter_scan_PMC_ks_kd(hostparams,nutrStar,paraRatio);
    
    if (q==1)
        plot(ppGpp_PMC(q,:),growthRate_PMC(q,:),'k-');
        errorbar(ZD2019_ppGpp_relA, ZD2019_Growth_relA, [], [], ZD2019_ppGpp_relA_err, ZD2019_ppGpp_relA_err, 'ko');
    end
    
    if (q==2)
        index1stzero = find(growthRate_PMC(q,:)==0, 1, 'first');
        plot(ppGpp_PMC(q,1:index1stzero-1),growthRate_PMC(q,1:index1stzero-1),'k-');
        plot(ppGpp_PMC(q,index1stzero-1:index1stzero),growthRate_PMC(q,index1stzero-1:index1stzero),'k--');
        plot(ppGpp_PMC(q,index1stzero:end),growthRate_PMC(q,index1stzero:end),'k-');
        errorbar(ZD2019_ppGpp_mesh, ZD2019_Growth_mesh, [], [], ZD2019_ppGpp_mesh_err, ZD2019_ppGpp_mesh_err, 'ko');
    end
    
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