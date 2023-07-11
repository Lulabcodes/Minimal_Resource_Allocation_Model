function run_ppGpp_dynamics(hp,nutr_pre,nutr_post)

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
                            'RelTol',tol,...
                            'AbsTol',tol,...
                            'Events',@myEvent_PMC);
%   time span                        
tspan = [0 10^10];

%--------------------
%   nutrient upshift
%--------------------

%   initial condition
x0 = [10,10,10];

[t1,x1] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr_pre(1),0,hp);
ppGpp1 = x1(:,3);
grate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,grate1(i)]  = Ecoli_GR_ODE_PMC(t1(end),x1(i,:),nutr_pre(1),0,hp);
end

[t2,x2] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x1(end,:),options_ode15s_PMC,nutr_post(1),0,hp);
ppGpp2 = x2(:,3);
grate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,grate2(i)]  = Ecoli_GR_ODE_PMC(t2(end),x2(i,:),nutr_post(1),0,hp);
end

figure();

%   maker size
MS = 8;

%   line width
LW = 1.5;

subplot(1,2,1);
hold on;

plot([t1-t1(end);t2]*60,[ppGpp1/ppGpp1(end);ppGpp2/ppGpp1(end)],'k-','LineWidth',LW);
axis square;
box on;
axis([-5,80,0,1]);
xlabel('Time (min)');
ylabel('ppGpp (\muM)');

%   Paper: Synthesis and turnover of basal level guanosine tetraphosphate
%   in Escherichia coli by J Friesen et al. (1975)
%   Shift from Tris-0.2% acetate (135 min) to rich medium (32 min)
timeUpshift     = [0,0.057,0.158,0.055,0.546,0.628,1.54,2.947,4.658,9.649,14.808,19.388,24.660,29.565,34.978,39.684,44.706,49.951,54.778,60.220,65.036,70.068]; %   min
ppGppUpshift    = [40.660,35.155,34.424,27.035,7.014,0.353,0.110,0.029,0.207,0.327,0.174,0.242,0.016,0.034,1.361,2.783,6.515,8.609,7.147,9.049,9.885,9.580];    %   pmol/A450
plot(timeUpshift,ppGppUpshift/ppGppUpshift(1),'ko','MarkerFaceColor',[61,191,255]/255,'MarkerSize',MS);

%---------------------
%   nutrient downshift
%---------------------

%   initial condition
x0 = [10,10,10];

[t1,x1] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutr_pre(2),0,hp);
ppGpp1 = x1(:,3);
grate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,grate1(i)]  = Ecoli_GR_ODE_PMC(t1(end),x1(i,:),nutr_pre(2),0,hp);
end

[t2,x2] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x1(end,:),options_ode15s_PMC,nutr_post(2),0,hp);
ppGpp2 = x2(:,3);
grate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,grate2(i)]  = Ecoli_GR_ODE_PMC(t2(end),x2(i,:),nutr_post(2),0,hp);
end

subplot(1,2,2);
hold on;

plot([t1-t1(end);t2]*60,[ppGpp1/ppGpp1(end);ppGpp2/ppGpp1(end)],'k-','LineWidth',LW);
axis square;
box on;
axis([-5,60,0,15]);
xlabel('Time (min)');
ylabel('ppGpp (\muM)');

%   Paper: Control of rRNA expression by small molecules is dynamic and
%   nonredundant by Murray et al. (2003)
%   
timeDownshift   = [0,1,3,5,10,20];    %   min
ppGppDownshift  = [1,2.90,8.71,12.45,10.16,5.02];    %   relative level
sdppGppDownshift= [0,0.33,0,0.78,0,0];
errorbar(timeDownshift,ppGppDownshift,sdppGppDownshift,'ko','MarkerSize',MS,'MarkerFaceColor',[61,191,255]/255);

end

