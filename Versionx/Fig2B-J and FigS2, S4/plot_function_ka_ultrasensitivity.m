function plot_function_ka_ultrasensitivity(kaScan,ppGppList,ppGppSynRate_kascan,growthRate_kascan,...
    kaList,res_PMC,res_TOC)

%   determine the index of ppGppList where growth rate is zero
k_nogrowth = length(kaScan);
for k=1:length(kaScan)
    if (pchip(kaList,res_PMC(:,4),kaScan(k)) <1e-3)
        k_nogrowth = k;
        break;
    end
end

%  Find index of ppGpp for maximum growth rate, lower and upper limits
IndexMaxGrowth_kascan  = zeros(1,length(kaScan));
IndexLowerLimit_kascan = zeros(1,length(kaScan));
IndexUpperLimit_kascan = zeros(1,length(kaScan));

for k=1:k_nogrowth
    
    [~,IndexMaxGrowth_kascan(k)]  = max(growthRate_kascan(k,:));
    try
        IndexLowerLimit_kascan(k)     = find(growthRate_kascan(k,1:IndexMaxGrowth_kascan(k)),1,'first');
        IndexUpperLimit_kascan(k)     = find(growthRate_kascan(k,IndexMaxGrowth_kascan(k):end),1,'last')+IndexMaxGrowth_kascan(k)-1;
    catch
        k_nogrowth = k;
        break;
    end
end

%   line width
LW  = 1.0;

%   color code
ColorPalette = [0,191,255;135,206,235;211,211,211;255,239,213]/255;

%-----------------------------
%   Growth rate vs. Ka
%-----------------------------

LS_maxgrowth    = zeros(1,length(kaScan));
for k=1:k_nogrowth-1

    ilow = IndexLowerLimit_kascan(k);
    iupp = IndexUpperLimit_kascan(k);
    deriv               = fnder(spline(ppGppList(ilow:iupp),ppGppSynRate_kascan(k,ilow:iupp)),1);
    localSensitivity    = abs(ppval(deriv,ppGppList(ilow:iupp)).*ppGppList(ilow:iupp)./ppGppSynRate_kascan(k,(ilow:iupp)));
    LS_maxgrowth(k)     = pchip(ppGppList(ilow:iupp),localSensitivity,ppGppList(IndexMaxGrowth_kascan(k)));
end

figure();
hold on;

growthRate_PMC_interp = pchip(kaList,res_PMC(:,4),kaScan(1:k_nogrowth-1));
growthRate_TOC_interp = pchip(kaList,res_TOC(:,4),kaScan(1:k_nogrowth-1));

plot(LS_maxgrowth(1:k_nogrowth-1),growthRate_PMC_interp,LS_maxgrowth(1:k_nogrowth-1),growthRate_TOC_interp);
axis square;
box on;
xlim([10^0,10^(2.5)]);
ylim([0,1.2]);
set(gca,'XTick',[10^0,10^(1.25),10^(2.5)]);
set(gca,'XTicklabel',{'10^{0}';'10^{1.25}';'10^{2.5}'});
set(gca,'YTick',[0,0.6,1.2]);
set(gca,'XScale','log');
xlabel('Sensitivity of ppGpp syn. rate at max growth');
ylabel('Growth rate (h^{-1})');

legend('PMC','TOC');
localS_WT= pchip(kaScan(1:k_nogrowth-1),LS_maxgrowth(1:k_nogrowth-1),10);
plot([localS_WT,localS_WT],[0,2],'k--');

figure();

%subplot(1,2,1);
hold on;

%plot(kaScan,maxLS,'-','Color',colorPalettes(1,:),'LineWidth',LW);
plot(kaScan(1:k_nogrowth-1),LS_maxgrowth(1:k_nogrowth-1),'-','Color',ColorPalette(2,:),'LineWidth',LW);
axis square;
box on;
xlabel('K_a (\muM)');
ylabel({'Local sensitivity'});
axis([min(kaScan),max(kaScan),0,100]);
set(gca,'XTick',[min(kaScan),10^((log10(max(kaScan))+log10(min(kaScan)))/2),max(kaScan)]);
set(gca,'YTick',[0,50,100]);
set(gca,'XScale','log');
legend('Sensitivity at maximal growth rate');

hpatch = patch([kaScan(k_nogrowth),kaScan(end),kaScan(end),kaScan(k_nogrowth)],...
    [0,0,100,100],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
uistack(hpatch,'bottom');

% %-----------------------------
% %   Growth rate vs. Ka
% %-----------------------------
% 
% subplot(1,2,2);
% hold on;
% 
% plot(kaList,res_PMC(:,4),'k-','Color',ColorPalette(1,:),'LineWidth',LW);
% plot(kaList,res_TOC(:,4),'k--','Color',ColorPalette(1,:),'LineWidth',LW);
% axis square;
% box on;
% xlabel('K_a (\muM)');
% ylabel('Growth rate (h^{-1})');
% 
% axis([1e-1,1e+7,0,1.2]);
% set(gca,'XScale','log');
% legend('PMC','TOC');
% 
% set(gca,'XTick',[1e-1,1e+3,1e+7]);
% set(gca,'YTick',[0,0.6,1.2]);

end