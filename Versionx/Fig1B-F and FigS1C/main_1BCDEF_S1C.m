clear;
clc;

addpath('../');

%%  read parameters
hostparams = readParameters();

%%  Simulation 1: Nutrient limitation
%   Fig1B
[nutr,growthRate] = run_nutrient_limitation(hostparams);
[growthRate,I] = unique(growthRate);
nutr           = nutr(I);

%%  Simulation 2: Translation inhibition
%   Fig1D
nutrList = [280,68,38,22];
run_translation_limitation(hostparams,nutrList);

%%  Simulation 3: Growth rate dynamics after nutrient upshift
%   Fig1E(left) and FigS1C
nutr_pre  = pchip(growthRate,nutr,[0.91,0.87,0.68,0.63,0.44,0.36,0.34,0.30,0.10]);
nutr_post = pchip(growthRate,nutr,[2.45,2.45,2.45,2.45,2.45,2.45,2.45,2.45,2.45]);
run_growthRate_dynamics(hostparams,nutr_pre,nutr_post);
 
%%  Simulation 4: ppGpp dynamics after nutrient upshift/downshift
%   Fig1E(right) and Fig1F
nutr_pre  = pchip(growthRate,nutr,[0.30,0.82]);
nutr_post = pchip(growthRate,nutr,[1.45,0.50]);
run_ppGpp_dynamics(hostparams,nutr_pre,nutr_post);

%%  Simulation 5: SpoT mutation
%   Fig1C

%   For succinate shift-up, we assume the ppGpp decay rate is measured for
%   the post-shift growth medium because the parameter came from Figure 4, 
%   where only ppGpp kinetic profiles after shift-up (t>0) are shown
kdPpGppList = [0.12,0.25]*60;   %   1/h
nutrList = pchip(growthRate,nutr,[1.28,1.50]*log(2));
run_SpoT_mutation(hostparams,kdPpGppList,nutrList);
