function run_SpoT_mutation(hp,kdPpGppList,nutrList)

%   hp: host cell parameters
%   kdPpGppList: ppGpp decay rate measured at different growth conditions
%   nutrList: nutrient level for the same growth conditions

growthRate = zeros(1,length(nutrList));

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3],...
                            'RelTol',tol,...
                            'AbsTol',tol,...
                            'Events',@myEvent_PMC);

%   time span                        
tspan = [0 10^10];

%   initial condition
x0 = [10,10,10];

for k=1:length(nutrList)
    
    hp_copy = hp;
    hp_copy.('kdPpGpp') = kdPpGppList(k);

    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,nutrList(k),0,hp_copy);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    [~,growthRate(k)] = Ecoli_GR_ODE_PMC(0,x(end,:),nutrList(k),0,hp_copy);
end

%   plot
figure();
hold on;

%   maker size
MS = 8;

%   bar width
barW = 0.4;

%   Growth rate of spoT mutants
%   Paper: Interaction of alleles of the relA, relC and spoT genes in
%   Escherichia coli: analysis of the interconversion of GTP, ppGpp, and
%   pppGpp
growthRate_spoT = [1.07,1.30]*log(2);

hb = bar(growthRate_spoT,barW,'FaceColor',[61,191,255]/255);
for ib=1:numel(hb.XData)
    plot(hb.XData(ib),growthRate(ib),'kx','MarkerSize',MS);
end
axis square;
box on;
axis([0,3,0,1]);
ylabel('Growth rate (h^{-1})');
set(gca,'YTick',[0,0.5,1]);
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{'Glucose';'Succinate'});
end

