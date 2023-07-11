function plot_cell_growth_vary_ppGpp(ppGppList,AminoAcid,Ribosome,growthRate,...
    ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth,hp)

%-----------------------------------------
%   Plot cell responses to ppGpp variation
%-----------------------------------------

figure();
hold on;

%   color code
ColorPalette = [0,191,255;135,206,235;211,211,211;255,239,213]/255;

%   line width
LW = 1.5;

for Flag = 0:0
    
    %   subplot(4,4,Flag+1);
    subplot(4,1,1);
    plot(ppGppList(IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),growthRate(Flag+1,IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),'LineWidth',LW,'Color',ColorPalette(1,:));
    hold on;
    
    if (IndexLowerLimit(Flag+1) > 1)
        plot(ppGppList(1:IndexLowerLimit(Flag+1)-1),growthRate(Flag+1,1:IndexLowerLimit(Flag+1)-1),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2,(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2],[0,growthRate(Flag+1,IndexLowerLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
    
    if (IndexUpperLimit(Flag+1) < length(ppGppList))
        plot(ppGppList(IndexUpperLimit(Flag+1)+1:end),growthRate(Flag+1,IndexUpperLimit(Flag+1)+1:end),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2,(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2],[0,growthRate(Flag+1,IndexUpperLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
    
    plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,1.6],'k--');
    axis square;
    box on;
    
    switch Flag
        case {0,1,2}
            axis([0,200,0,1.2]);
            set(gca,'XTick',[0,100,200]);
        otherwise
            axis([0,600,0,1.2]);
            set(gca,'XTick',[0,300,600]);
    end
    
    set(gca,'YTick',[0,0.6,1.2]);
    ylabel('Growth rate (h^{-1})');
    hl = patch([ppGppList(1),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(1)],...
        [0,0,2.0,2.0],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1))],...
        [0,0,2.0,2.0],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList(IndexUpperLimit(Flag+1)),ppGppList(end),ppGppList(end),ppGppList(IndexUpperLimit(Flag+1))],...
        [0,0,2.0,2.0],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
%     switch Flag
%         case 0
%             title('Both saturated');
%         case 1
%             title('tRNA charging saturated');
%         case 2
%             title('Translation saturated');
%         case 3
%             title('Both non-saturated');
%     end
    
    %   subplot(4,4,Flag+5);
    subplot(4,1,2);
    plot(ppGppList(IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),Ribosome(Flag+1,IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),'LineWidth',LW,'Color',ColorPalette(1,:));
    hold on;
    
    if (IndexLowerLimit(Flag+1) > 1)
        plot(ppGppList(1:IndexLowerLimit(Flag+1)-1),Ribosome(Flag+1,1:IndexLowerLimit(Flag+1)-1),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2,(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2],[Ribosome(Flag+1,IndexLowerLimit(Flag+1)-1),Ribosome(Flag+1,IndexLowerLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
    
    if (IndexUpperLimit(Flag+1) < length(ppGppList))
        plot(ppGppList(IndexUpperLimit(Flag+1)+1:end),Ribosome(Flag+1,IndexUpperLimit(Flag+1)+1:end),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2,(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2],[Ribosome(Flag+1,IndexUpperLimit(Flag+1)+1),Ribosome(Flag+1,IndexUpperLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
    
    plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,300],'k--');
    axis square;
    box on;
    
    switch Flag
        case {0,1,2}
            axis([0,200,0,150]);
            set(gca,'XTick',[0,100,200]);
        otherwise
            axis([0,600,0,150]);
            set(gca,'XTick',[0,300,600]);
    end
    
    set(gca,'YTick',[0,75,150]);
    ylabel('Ribosome (\muM)');
    hl = patch([ppGppList(1),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(1)],...
        [0,0,300.0,300.0],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1))],...
        [0,0,300.0,300.0],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList(IndexUpperLimit(Flag+1)),ppGppList(end),ppGppList(end),ppGppList(IndexUpperLimit(Flag+1))],...
        [0,0,300.0,300.0],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
    
    %   subplot(4,4,Flag+9);
    subplot(4,1,3);
    plot(ppGppList(IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),AminoAcid(Flag+1,IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),'LineWidth',LW,'Color',ColorPalette(1,:));
    hold on;
    
%     if (IndexLowerLimit(Flag+1) > 1)
%         plot(ppGppList(1:IndexLowerLimit(Flag+1)-1),AminoAcid(Flag+1,1:IndexLowerLimit(Flag+1)-1),'LineWidth',LW,'Color',ColorPalette(1,:));
%         plot([(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2,(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2],[AminoAcid(Flag+1,IndexLowerLimit(Flag+1)-1),AminoAcid(Flag+1,IndexLowerLimit(Flag+1))],'--','Color',ColorPalette(1,:));
%     end
%     
%     if (IndexUpperLimit(Flag+1) < length(ppGppList))
%         plot(ppGppList(IndexUpperLimit(Flag+1)+1:end),AminoAcid(Flag+1,IndexUpperLimit(Flag+1)+1:end),'LineWidth',LW,'Color',ColorPalette(1,:));
%         plot([(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2,(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2],[AminoAcid(Flag+1,IndexUpperLimit(Flag+1)+1),AminoAcid(Flag+1,IndexUpperLimit(Flag+1))],'--','Color',ColorPalette(1,:));
%     end
    
    plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[1e-1,1e+5],'k--');
    axis square;
    box on;
    
    switch Flag
        case {0,1,2}
            axis([0,200,1e-1,1e+5]);
            set(gca,'XTick',[0,100,200]);
        otherwise
            axis([0,600,1e-1,1e+5]);
            set(gca,'XTick',[0,300,600]);
    end
    
    set(gca,'YScale','log');
    set(gca,'YTick',10.^[-1,2,5]);
    set(gca,'YTickLabel',{'10^{-1}';'10^{2}';'10^{5}'});
    ylabel('Amino acids (\muM)');
    hl = patch([ppGppList(1),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(1)],...
        [1e-1,1e-1,1e+5,1e+5],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1))],...
        [1e-1,1e-1,1e+5,1e+5],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList(IndexUpperLimit(Flag+1)),ppGppList(end),ppGppList(end),ppGppList(IndexUpperLimit(Flag+1))],...
        [1e-1,1e-1,1e+5,1e+5],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
    
    %   subplot(4,4,Flag+13);
    subplot(4,1,4);
    plot(ppGppList(IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1):IndexUpperLimit(Flag+1)),'LineWidth',LW,'Color',ColorPalette(1,:));
    hold on;

    if (IndexLowerLimit(Flag+1) > 1)
        plot(ppGppList(1:IndexLowerLimit(Flag+1)-1),ppGppSynRate(Flag+1,1:IndexLowerLimit(Flag+1)-1),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2,(ppGppList(IndexLowerLimit(Flag+1)-1)+ppGppList(IndexLowerLimit(Flag+1)))/2],[ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1)-1),ppGppSynRate(Flag+1,IndexLowerLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
        
    if (IndexUpperLimit(Flag+1) < length(ppGppList))
        plot(ppGppList(IndexUpperLimit(Flag+1)+1:end),ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1)+1:end),'LineWidth',LW,'Color',ColorPalette(1,:));
        plot([(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2,(ppGppList(IndexUpperLimit(Flag+1))+ppGppList(IndexUpperLimit(Flag+1)+1))/2],[ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1)+1),ppGppSynRate(Flag+1,IndexUpperLimit(Flag+1))],'--','Color',ColorPalette(1,:));
    end
    
    %    plot(ppGppList,ppGppDegRate(Flag+1,:),'--','LineWidth',LW,'Color',ColorPalette(1,:));
    %legend({'syn.','deg.'},'Location','SouthWest');
    plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[1e1,1e+7],'k--');
    axis square;
    box on;
    
    switch Flag
        case {0,1,2}
            axis([0,200,1e+1,1e+7]);
            set(gca,'XTick',[0,100,200]);
        otherwise
            axis([0,600,1e+1,1e+7]);
            set(gca,'XTick',[0,300,600]);
    end
    
    set(gca,'YScale','log');
    set(gca,'YTick',10.^[1,4,7]);
    set(gca,'YTickLabel',{'10^{1}';'10^{4}';'10^{7}'});
    xlabel('ppGpp (\muM)');
    ylabel('dppGpp/dt (\muM/h)');
    hl = patch([ppGppList(1),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1)),ppGppList(1)],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    hm = patch([ppGppList(IndexLowerLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexUpperLimit(Flag+1)),ppGppList(IndexLowerLimit(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(4,:),'EdgeColor','none','FaceAlpha',0.5);
    hr = patch([ppGppList(IndexUpperLimit(Flag+1)),ppGppList(end),ppGppList(end),ppGppList(IndexUpperLimit(Flag+1))],...
        [1e-2,1e-2,1e+7,1e+7],ColorPalette(3,:),'EdgeColor','none','FaceAlpha',0.5);
    uistack(hl,'bottom');
    uistack(hm,'bottom');
    uistack(hr,'bottom');
end

end

