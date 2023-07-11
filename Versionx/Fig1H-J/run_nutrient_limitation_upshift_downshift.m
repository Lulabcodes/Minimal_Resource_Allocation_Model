function [t_PMC_upshift,growthRate_PMC_t_upshift,ppGpp_PMC_t_upshift,...
    t_PMC_downshift,growthRate_PMC_t_downshift,ppGpp_PMC_t_downshift,...
    t_TOC_upshift,growthRate_TOC_t_upshift,ppGpp_TOC_t_upshift,...
    t_TOC_downshift,growthRate_TOC_t_downshift,ppGpp_TOC_t_downshift,...
    t_CC_upshift,growthRate_CC_t_upshift,ppGpp_CC_t_upshift,...
    t_CC_downshift,growthRate_CC_t_downshift,ppGpp_CC_t_downshift] = ...
    run_nutrient_limitation_upshift_downshift(hp,shiftNutr,shiftPpGpp_TOC,shiftPpGpp_CC)

%   hp: host cell parameters
%   nutrList: list of nutrient values

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC);
options_ode15s_CC = odeset('NonNegative',[1,2],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_CC);

%   time span
tspan = [0 10^10];

%%  ppGpp-mediated control

%   upshift
x0 = [10,10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,shiftNutr(3),0,hp);
ppGpp1      = x1(:,3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_PMC(t1(end),x1(i,:),shiftNutr(3),0,hp);
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x1(end,:),options_ode15s_PMC,shiftNutr(2),0,hp);
ppGpp2      = x2(:,3);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_PMC(t2(end),x2(i,:),shiftNutr(2),0,hp);
end
t_PMC_upshift = [t1-t1(end);t2];
ppGpp_PMC_t_upshift = [ppGpp1;ppGpp2];
growthRate_PMC_t_upshift = [growthRate1';growthRate2'];

%   downshift
x0 = [10,10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,shiftNutr(3),0,hp);
ppGpp1      = x1(:,3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_PMC(t1(end),x1(i,:),shiftNutr(3),0,hp);
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x1(end,:),options_ode15s_PMC,shiftNutr(1),0,hp);
ppGpp2      = x2(:,3);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_PMC(t2(end),x2(i,:),shiftNutr(1),0,hp);
end
t_PMC_downshift = [t1-t1(end);t2];
ppGpp_PMC_t_downshift = [ppGpp1;ppGpp2];
growthRate_PMC_t_downshift = [growthRate1';growthRate2'];


%%  theoretical optimal control

%   upshift
x0 = [10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,shiftNutr(3),0,hp,shiftPpGpp_TOC(3));
ppGpp1      = ones(length(t1),1)*shiftPpGpp_TOC(3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_CC(t1(end),x1(i,:),shiftNutr(3),0,hp,shiftPpGpp_TOC(3));
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_CC,tspan,x1(end,:),options_ode15s_CC,shiftNutr(2),0,hp,shiftPpGpp_TOC(2));
ppGpp2      = ones(length(t2),1)*shiftPpGpp_TOC(2);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_CC(t2(end),x2(i,:),shiftNutr(2),0,hp,shiftPpGpp_TOC(2));
end
t_TOC_upshift = [t1-t1(end);t2];
ppGpp_TOC_t_upshift = [ppGpp1;ppGpp2];
growthRate_TOC_t_upshift = [growthRate1';growthRate2'];

%   downshift
x0 = [10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,shiftNutr(3),0,hp,shiftPpGpp_TOC(3));
ppGpp1      = ones(length(t1),1)*shiftPpGpp_TOC(3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_CC(t1(end),x1(i,:),shiftNutr(3),0,hp,shiftPpGpp_TOC(3));
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_CC,tspan,x1(end,:),options_ode15s_CC,shiftNutr(1),0,hp,shiftPpGpp_TOC(1));
ppGpp2      = ones(length(t2),1)*shiftPpGpp_TOC(1);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_CC(t2(end),x2(i,:),shiftNutr(1),0,hp,shiftPpGpp_TOC(1));
end
t_TOC_downshift = [t1-t1(end);t2];
ppGpp_TOC_t_downshift = [ppGpp1;ppGpp2];
growthRate_TOC_t_downshift = [growthRate1';growthRate2'];

%%  constant control

%   upshift
x0 = [10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,shiftNutr(3),0,hp,shiftPpGpp_CC(3));
ppGpp1      = ones(length(t1),1)*shiftPpGpp_CC(3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_CC(t1(end),x1(i,:),shiftNutr(3),0,hp,shiftPpGpp_CC(3));
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_CC,tspan,x1(end,:),options_ode15s_CC,shiftNutr(2),0,hp,shiftPpGpp_CC(2));
ppGpp2      = ones(length(t2),1)*shiftPpGpp_CC(2);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_CC(t2(end),x2(i,:),shiftNutr(2),0,hp,shiftPpGpp_CC(2));
end
t_CC_upshift = [t1-t1(end);t2];
ppGpp_CC_t_upshift = [ppGpp1;ppGpp2];
growthRate_CC_t_upshift = [growthRate1';growthRate2'];

%   downshift
x0 = [10,10];
[t1,x1] = ode15s(@Ecoli_GR_ODE_CC,tspan,x0,options_ode15s_CC,shiftNutr(3),0,hp,shiftPpGpp_CC(3));
ppGpp1      = ones(length(t1),1)*shiftPpGpp_CC(3);
growthRate1 = zeros(1,length(t1));
for i=1:length(t1)
    [~,growthRate1(i)]  = Ecoli_GR_ODE_CC(t1(end),x1(i,:),shiftNutr(3),0,hp,shiftPpGpp_CC(3));
end
[t2,x2] = ode15s(@Ecoli_GR_ODE_CC,tspan,x1(end,:),options_ode15s_CC,shiftNutr(1),0,hp,shiftPpGpp_CC(1));
ppGpp2      = ones(length(t2),1)*shiftPpGpp_CC(1);
growthRate2 = zeros(1,length(t2));
for i=1:length(t2)
    [~,growthRate2(i)]  = Ecoli_GR_ODE_CC(t2(end),x2(i,:),shiftNutr(1),0,hp,shiftPpGpp_CC(1));
end
t_CC_downshift = [t1-t1(end);t2];
ppGpp_CC_t_downshift = [ppGpp1;ppGpp2];
growthRate_CC_t_downshift = [growthRate1';growthRate2'];

end