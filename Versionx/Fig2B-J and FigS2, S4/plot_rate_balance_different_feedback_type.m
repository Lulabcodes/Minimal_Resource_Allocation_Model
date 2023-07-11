function plot_rate_balance_different_feedback_type(...
    ppGppList_fb,ppGppSynRate_fb,ppGppDegRate_fb,growthRate_fb,IndexLowerLimit_fb,IndexUpperLimit_fb)

%   line width
LW = 1.5;

%   color code
ColorPalette = [0,191,255;135,206,235;211,211,211;255,239,213]/255;


for Flag=0:2

    subplot(2,3,Flag+1);
    hold on;

    plot(ppGppList_fb(IndexLowerLimit_fb(Flag+1):IndexUpperLimit_fb(Flag+1)),ppGppSynRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1):IndexUpperLimit_fb(Flag+1)),'k-','LineWidth',LW);
    plot(ppGppList_fb,ppGppDegRate_fb(Flag+1,:),'k--');
    axis square;
    box on;
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    axis([1e-3,1e+5,1e-1,1e+7]);
    ylabel('dppGpp/dt (\muM/h)');
    xlabel('ppGpp (\muM)');
    
    if (IndexLowerLimit_fb(Flag+1) > 1)
        plot(ppGppList_fb(1:IndexLowerLimit_fb(Flag+1)-1),ppGppSynRate_fb(Flag+1,1:IndexLowerLimit_fb(Flag+1)-1),'k-','LineWidth',LW);
        plot([(ppGppList_fb(IndexLowerLimit_fb(Flag+1)-1)+ppGppList_fb(IndexLowerLimit_fb(Flag+1)))/2,(ppGppList_fb(IndexLowerLimit_fb(Flag+1)-1)+ppGppList_fb(IndexLowerLimit_fb(Flag+1)))/2],[ppGppSynRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1)-1),ppGppSynRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1))],'k--');
    end
        
    if (IndexUpperLimit_fb(Flag+1) < length(ppGppList_fb))
        plot(ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1:end),ppGppSynRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1)+1:end),'k','LineWidth',LW);
        plot([(ppGppList_fb(IndexUpperLimit_fb(Flag+1))+ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1))/2,(ppGppList_fb(IndexUpperLimit_fb(Flag+1))+ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1))/2],[ppGppSynRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1)+1),ppGppSynRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1))],'k--');
    end
    
    hl = patch([ppGppList_fb(1),ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(1)],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(IndexLowerLimit_fb(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(end),ppGppList_fb(end),ppGppList_fb(IndexUpperLimit_fb(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
    
    subplot(2,3,Flag+4);
    hold on;

    plot(ppGppList_fb(IndexLowerLimit_fb(Flag+1):IndexUpperLimit_fb(Flag+1)),growthRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1):IndexUpperLimit_fb(Flag+1)),'k-','LineWidth',LW);
    axis square;
    box on;
    set(gca,'XScale','log');
    axis([1e-3,1e+5,0,1.6]);
    ylabel('Growth rate (h^{-1})');
    xlabel('ppGpp (\muM)');
    
    if (IndexLowerLimit_fb(Flag+1) > 1)
        plot(ppGppList_fb(1:IndexLowerLimit_fb(Flag+1)-1),growthRate_fb(Flag+1,1:IndexLowerLimit_fb(Flag+1)-1),'k-','LineWidth',LW);
        plot([(ppGppList_fb(IndexLowerLimit_fb(Flag+1)-1)+ppGppList_fb(IndexLowerLimit_fb(Flag+1)))/2,(ppGppList_fb(IndexLowerLimit_fb(Flag+1)-1)+ppGppList_fb(IndexLowerLimit_fb(Flag+1)))/2],[growthRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1)-1),growthRate_fb(Flag+1,IndexLowerLimit_fb(Flag+1))],'k--');
    end
        
    if (IndexUpperLimit_fb(Flag+1) < length(ppGppList_fb))
        plot(ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1:end),ppGppSynRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1)+1:end),'k','LineWidth',LW);
        plot([(ppGppList_fb(IndexUpperLimit_fb(Flag+1))+ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1))/2,(ppGppList_fb(IndexUpperLimit_fb(Flag+1))+ppGppList_fb(IndexUpperLimit_fb(Flag+1)+1))/2],[growthRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1)+1),growthRate_fb(Flag+1,IndexUpperLimit_fb(Flag+1))],'k--');
    end
    
    hl = patch([ppGppList_fb(1),ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(1)],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList_fb(IndexLowerLimit_fb(Flag+1)),ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(IndexLowerLimit_fb(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList_fb(IndexUpperLimit_fb(Flag+1)),ppGppList_fb(end),ppGppList_fb(end),ppGppList_fb(IndexUpperLimit_fb(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
end

maxGrowth= max(growthRate_fb(1,:));

subplot(2,3,1);
hold on;
xlim([1e1,1e3]);
set(gca,'XTick',[1e+1,1e+2,1e+3]);
set(gca,'YTick',[1e-1,1e+3,1e+7]);

subplot(2,3,2);
xlim([1e2,1e4]);
set(gca,'XTick',[1e+2,1e+3,1e+4]);
set(gca,'YTick',[1e-1,1e+3,1e+7]);

subplot(2,3,3);
set(gca,'XTick',[1e-3,1e+1,1e+5]);
set(gca,'YTick',[1e-1,1e+3,1e+7]);

subplot(2,3,4);
hold on;
xlim([1e1,1e3]);
ylim([0,1.2]);
set(gca,'XTick',[1e+1,1e+2,1e+3]);
set(gca,'YTick',[0 0.6 1.2]);
plot([1e1,1e3],[maxGrowth,maxGrowth],'r--');

subplot(2,3,5);
hold on;
xlim([1e2,1e4]);
ylim([0,1.2]);
set(gca,'XTick',[1e+2,1e+3,1e+4]);
set(gca,'YTick',[0 0.6 1.2]);
plot([1e2,1e4],[maxGrowth,maxGrowth],'r--');

subplot(2,3,6);
hold on;
ylim([0,1.2]);
set(gca,'XTick',[1e-3,1e+1,1e+5]);
set(gca,'YTick',[0 0.6 1.2]);
plot([1e-3,1e5],[maxGrowth,maxGrowth],'r--');

end

