function plot_growth_rate_optimization_vary_Ka(kaList,res_PMC,res_TOC)

%-----------------------------------------
%   Plot cell responses to ppGpp variation
%-----------------------------------------

figure();

%   color code
ColorPalette = [255,69,0;34,139,34]/255;

%   line width
LW = 1.5;

subplot(1,2,1);
hold on;

plot(kaList,res_PMC(:,4),'k-','Color',ColorPalette(1,:),'LineWidth',LW);
plot(kaList,res_TOC(:,4),'k-','Color',ColorPalette(2,:),'LineWidth',LW);
axis square;
box on;
xlabel('K_a (\muM)');
ylabel('Growth rate (h^{-1})');

axis([min(kaList),max(kaList),0,1.2]);
set(gca,'XScale','log');
legend('PMC','TOC');

set(gca,'YTick',[0,0.6,1.2]);
set(gca,'XTick',[1e-1,1e+3,1e+7]);

subplot(1,2,2);
hold on;

plot(kaList,res_PMC(:,3),'k-','Color',ColorPalette(1,:),'LineWidth',LW);
plot(kaList,res_TOC(:,3),'k-','Color',ColorPalette(2,:),'LineWidth',LW);
axis square;
box on;
xlabel('K_a (\muM)');
ylabel('ppGpp (\muM)');

axis([min(kaList),max(kaList),0,150]);
set(gca,'XScale','log');
legend('PMC','TOC');

end

