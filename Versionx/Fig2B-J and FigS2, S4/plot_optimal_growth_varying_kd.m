function plot_optimal_growth_varying_kd(hp,ppGppList,growthRate,ppGppSynRate,IndexMaxGrowth,IndexLowerLimit,IndexUpperLimit)

%-------------------------------------------------------------------------
%   Plot ratio of growth rate to optimal growth rate for different models
%-------------------------------------------------------------------------

figure();
hold on;

%   line width
LW = 1.5;

%   Marker size
MS = 8;

%   ppGpp range for data processing
window = [10,1000;10,1000;10,1000;10,2000];

for Flag = 0:0
    
    [~,IndexL] = min(abs(ppGppList-window(Flag+1,1)));
    [~,IndexR] = min(abs(ppGppList-window(Flag+1,2)));

    %   subplot(1,4,Flag+1);
    hold on;
    
    %   PMC
    plot(ppGppList(IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),'LineWidth',LW,'Color',[0,191,255]/255);
    hold on;

    if (IndexLowerLimit(Flag+1) > 1)
        plot(ppGppList(1:IndexLowerLimit(Flag+1)-1),ppGppSynRate(Flag+1,1:IndexLowerLimit(Flag+1)-1),'LineWidth',LW,'Color',[0,191,255]/255);
        plot([(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2,(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2],[ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1)-1),ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1))],'--','Color',[0,191,255]/255);
    end
        
    if (IndexUpperLimit(Flag+1) < length(ppGppList))
        plot(ppGppList(IndexUpperLimit(Flag+1)+1:end),ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1)+1:end),'LineWidth',LW,'Color',[0,191,255]/255);
        plot([(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2,(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2],[ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1)+1),ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1))],'--','Color',[0,191,255]/255);
    end
    ppGppDegRate1    = (hp.('kdPpGpp')*0.1+growthRate(Flag+1,:)).*ppGppList;
    ppGppDegRate2    = (hp.('kdPpGpp')*1+growthRate(Flag+1,:)).*ppGppList;
    ppGppDegRate3    = (hp.('kdPpGpp')*10+growthRate(Flag+1,:)).*ppGppList;
    [~,IndexIntersec1]  = min(abs(ppGppSynRate(Flag+1,IndexL:IndexR)-ppGppDegRate1(IndexL:IndexR)));
    IndexIntersec1 = IndexIntersec1+IndexL-1;
    [~,IndexIntersec10] = min(abs(ppGppSynRate(Flag+1,IndexL:IndexR)-ppGppDegRate2(IndexL:IndexR)));
    IndexIntersec10 = IndexIntersec10+IndexL-1;
    [~,IndexIntersec100]= min(abs(ppGppSynRate(Flag+1,IndexL:IndexR)-ppGppDegRate3(IndexL:IndexR)));
    IndexIntersec100 = IndexIntersec100+IndexL-1;
    
    plot(ppGppList,ppGppDegRate1,'--','LineWidth',LW,'Color',[0,191,255]/255);
    plot(ppGppList,ppGppDegRate2,'--','LineWidth',LW,'Color',[0,191,255]/255);
    plot(ppGppList,ppGppDegRate3,'--','LineWidth',LW,'Color',[0,191,255]/255);
    plot(ppGppList(IndexIntersec1),ppGppSynRate(Flag+1,IndexIntersec1),'ko','MarkerSize',MS,'MarkerFaceColor',[0,191,255]/255);
    plot(ppGppList(IndexIntersec10),ppGppSynRate(Flag+1,IndexIntersec10),'ko','MarkerSize',MS,'MarkerFaceColor',[0,191,255]/255);
    plot(ppGppList(IndexIntersec100),ppGppSynRate(Flag+1,IndexIntersec100),'ko','MarkerSize',MS,'MarkerFaceColor',[0,191,255]/255);
    
    %   shading
    hl = patch([ppGppList(IndexIntersec100),ppGppList(IndexIntersec1),ppGppList(IndexIntersec1),ppGppList(IndexIntersec100)],...
        [1e1,1e1,1e+7,1e+7],[0,191,255]/255,'EdgeColor','none','FaceAlpha',0.3);
    uistack(hl,'bottom');
    
    %   PFC
    ppGppSynRate_PFC = ppGppSynRate(Flag+1,IndexIntersec10)*ppGppList(IndexIntersec10)./ppGppList;
    plot(ppGppList,ppGppSynRate_PFC,'LineWidth',LW,'Color',[128,128,128]/255);
    [~,IndexIntersec1]  = min(abs(ppGppSynRate_PFC-ppGppDegRate1));
    [~,IndexIntersec100]= min(abs(ppGppSynRate_PFC-ppGppDegRate3));
    plot(ppGppList(IndexIntersec1),ppGppSynRate_PFC(IndexIntersec1),'ko','MarkerSize',MS,'MarkerFaceColor',[128,128,128]/255);
    plot(ppGppList(IndexIntersec100),ppGppSynRate_PFC(IndexIntersec100),'ko','MarkerSize',MS,'MarkerFaceColor',[128,128,128]/255);
    
    %   shading
    hl = patch([ppGppList(IndexIntersec100),ppGppList(IndexIntersec1),ppGppList(IndexIntersec1),ppGppList(IndexIntersec100)],...
        [1e1,1e1,1e+7,1e+7],[128,128,128]/255,'EdgeColor','none','FaceAlpha',0.3);
    uistack(hl,'bottom');
    
    %legend({'syn.','deg.'},'Location','SouthWest');
    plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[1e+1,1e+8],'k--');
    axis square;
    box on;
    set(gca,'YScale','log');
    axis([10,1000,1e+1,1e+7]);
    set(gca,'XScale','log');
    set(gca,'XTick',10.^[0:1:3]);
    set(gca,'XTickLabel',{'10^{0}';'10^{1}';'10^{2}';'10^{3}'});
    set(gca,'YTick',10.^[1:3:7]);
    set(gca,'YTickLabel',{'10^{1}';'10^{4}';'10^{7}'});
    xlabel('ppGpp (\muM)');
    ylabel('ppGpp syn./deg. rate (\muM/h)');
    
end

end

