function run_growthRate_dynamics(hp,nutr_pre,nutr_post)

%   Model by Mori et al 2017
%   Mori_lambda0 = [1.04,0.99,0.86,0.77,0.72,0.64,0.68,0.66,0.5];
%   growthRateMori2017 = Mori_lambda0(k)+(2.45-Mori_lambda0(k))./(1+(2.45/Mori_lambda0(k))./(exp(2.45*t2)-1));
%   growthRateMori2017 = growthRateMori2017';
    
%   The following data were extracted from the same paper
%   Paper: Quantifying the benefit of a proteome reserve in fluctuating
%   environments by Matteo Mori 2017

%   0.2% Glucose
Time_Glu02      = [-16.42,7.36,11.89,16.98,22.08,27.17,31.98,37.08,42.17,46.98,52.08];
GR_Glu02        = [0.96,0.95,1.29,1.54,1.75,1.93,2.01,1.93,1.86,1.98,1.96];
GR_Glu02_err    = [0.09,0.19,0.19,0.19,0.20,0.20,0.19,0.19,0.19,0.20,0.18];

%   0.2% Arabinose
Time_Arab02     = [-25.03,-7.39,6.88,12.90,18.14,22.66,27.83,32.80,37.80,43.34,48.51,53.19];
GR_Arab02       = [0.86,0.84,0.87,1.66,1.76,1.48,1.51,1.37,1.51,1.87,1.92,1.77];
GR_Arab02_err   = [0.11,0.10,0.29,0.28,0.28,0.29,0.30,0.29,0.29,0.28,0.28,0.29];

%   0.2% Glycerol
Time_Gly02      = [-14.64,7.92,12.78,18.17,23.51,28.27,33.46,38.19,43.25,48.63,53.76];
GR_Gly02        = [0.72,0.94,0.99,1.32,1.59,1.54,1.64,1.55,1.51,1.83,1.87];
GR_Gly02_err    = [0.11,0.14,0.15,0.15,0.14,0.15,0.15,0.16,0.14,0.15,0.15];

%   0.2% Fructose
Time_Fruc02     = [-9.05,5.59,10.27,15.19,19.74,24.71,29.77,34.60,39.60,44.39,49.25];
GR_Fruc02       = [0.60,0.69,0.84,1.02,1.32,1.45,1.50,1.47,1.58,1.91,1.86];
GR_Fruc02_err   = [0.12,0.15,0.15,0.17,0.15,0.16,0.16,0.15,0.15,0.16,0.15];

%   0.2% Mannose
Time_Mann02     = [-25.42,4.19,8.86,14.02,17.85,22.95,27.04,31.57,36.70,41.36,46.66];
GR_Mann02       = [0.46,0.70,0.84,0.82,1.14,1.14,1.46,1.56,1.55,1.69,1.54];
GR_Mann02_err   = [0.10,0.14,0.15,0.16,0.13,0.14,0.15,0.14,0.14,0.14,0.14];

%   0.1% Mannose
Time_Mann01     = [-15.26,5.55,10.50,15.16,19.93,24.85,30.01,34.92,40.11,45.48];
GR_Mann01       = [0.38,0.57,0.78,1.01,1.48,1.32,1.23,1.49,1.77,1.79];
GR_Mann01_err   = [0.10,0.16,0.16,0.16,0.16,0.17,0.15,0.15,0.14,0.16];

%   0.075% Mannose
Time_Mann0075   = [-16.46,9.70,14.33,19.22,25.96,32.81,37.67,41.36,47.64,53.67,58.07];
GR_Mann0075     = [0.30,0.75,0.70,0.64,1.12,1.64,1.57,1.24,1.59,1.95,1.74];
GR_Mann0075_err = [0.10,0.20,0.20,0.20,0.21,0.20,0.20,0.20,0.20,0.19,0.20];

%   0.05% Mannose
Time_Mann005    = [-17.79,6.47,11.32,16.36,21.12,25.34,30.87,35.84,40.47,45.28,50.38];
GR_Mann005      = [0.26,0.66,0.93,0.68,1.05,1.70,1.22,1.36,1.88,1.87,1.87];
GR_Mann005_err  = [0.10,0.22,0.24,0.23,0.23,0.22,0.22,0.21,0.22,0.23,0.23];

%   20mM Glutamate
Time_Glut20     = [7.53,12.89,18.17,23.39,29.04,34.03,38.78,44.12,49.58,54.87];
GR_Glut20       = [0.34,0.69,0.80,0.88,1.23,1.35,1.13,1.28,1.69,1.82];
GR_Glut20_err   = [0.10,0.10,0.11,0.11,0.09,0.10,0.10,0.10,0.11,0.10];

%   initial condition
x0 = [10,10,10];

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
                            'RelTol',tol,...
                            'AbsTol',tol,...
                            'Events',@myEvent_PMC);
%   time span                        
tspan = [0 10^10];

figure();

%   maker size
MS = 8;

%   line width
LW = 1.5;

for k=1:length(nutr_pre)

    [t1,x1] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr_pre(k),0,hp);
    growthRate1 = zeros(1,length(t1));
    for i=1:length(t1)
        [~,growthRate1(i)]  = Ecoli_GR_ODE_PMC(t1(end),x1(i,:),nutr_pre(k),0,hp);
    end
    
    [t2,x2] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x1(end,:),options_ode15s_PMC,nutr_post(k),0,hp);
    growthRate2 = zeros(1,length(t2));
    for i=1:length(t2)
        [~,growthRate2(i)]  = Ecoli_GR_ODE_PMC(t2(end),x2(i,:),nutr_post(k),0,hp);
    end

    switch k
        case 1
            subplot(3,3,1);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Glu02,GR_Glu02,GR_Glu02_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,1.04,0.10,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
            
        case 2
            subplot(3,3,2);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Arab02,GR_Arab02,GR_Arab02_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.99,0.10,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 3
            subplot(3,3,3);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Gly02,GR_Gly02,GR_Gly02_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.86,0.09,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 4
            subplot(3,3,4);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Fruc02,GR_Fruc02,GR_Fruc02_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.77,0.08,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 5
            subplot(3,3,5);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Mann02,GR_Mann02,GR_Mann02_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.72,0.08,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 6
            subplot(3,3,6);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Mann01,GR_Mann01,GR_Mann01_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.64,0.08,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 7
            subplot(3,3,7);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Mann0075,GR_Mann0075,GR_Mann0075_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.68,0.08,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 8
            subplot(3,3,8);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Mann005,GR_Mann005,GR_Mann005_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.66,0.07,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
        case 9
            subplot(3,3,9);
            plot([t1-t1(end);t2]*60,[growthRate1';growthRate2'],'k-','LineWidth',LW);
            hold on;
            %   plot([t1-t1(end);t2]*60,[growthRate1';growthRateMori2017'],'r-');
            errorbar(Time_Glut20,GR_Glut20,GR_Glut20_err,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);
            axis square;
            box on;
            axis([-30,60,0,2.5]);
            %errorbar(0,0.5,0.06,'ro','MarkerSize',8);
            xlabel('Time (h)');
            ylabel('Growth rate (h^{-1})');
    end
%     [nutr_post(k)*7336*1.6/325/19/3600*1.1,0.05]
%     plot([0],[growthRate2(end)*x1(end,2)/x2(end,2)*(x2(end,2)*0.5-x1(end,2))/(x2(end,2)-x1(end,2))],'x');
end



end

