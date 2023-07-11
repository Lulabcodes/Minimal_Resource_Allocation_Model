function plot_rate_balance_analysis_full_model_vary_Ka(...
    ppGppToShow,AminoAcidToShow,AminoAcidSynRate,AminoAcidDegRate,ppGppList,AminoAcid,IndexMaxGrowth)

figure();

%   line width
LW  = 1.0;

%   marker size
MS  = 6;

%   color code
cc = jet(length(ppGppToShow));

Flag=0;

%   plot synthesis and degradation curves

subplot(2,1,1);
hold on;

AminoAcidIntersect = zeros(1,length(ppGppToShow));
for k=1:length(ppGppToShow)
    
    plot(AminoAcidToShow,squeeze(AminoAcidSynRate(Flag+1,k,:)),'k-','linewidth',LW,'Color',cc(k,:));
    plot(AminoAcidToShow,squeeze(AminoAcidDegRate(Flag+1,k,:)),'k--','linewidth',LW,'Color',cc(k,:));
    
    %    x and y coordinates of intersection of amino acid
    %    synthesis/degradation rates
    xInterp    = pchip(ppGppList,AminoAcid(Flag+1,:),ppGppToShow(k));
    yInterp    = pchip(AminoAcidToShow,squeeze(AminoAcidSynRate(Flag+1,k,:)),xInterp);
    AminoAcidIntersect(k) = xInterp;
    plot(xInterp,yInterp,'ko','MarkerSize',MS,'MarkerFaceColor',cc(k,:));
end
axis square;
box on;
if (Flag == 3)
    axis([1e+0,1e+4,0,3e+5]);
    set(gca,'XTick',[1e+0,1e+2,1e+4]);
    set(gca,'XTicklabel',{'10^{0}';'10^{2}';'10^{4}'});
else
    axis([1e-1,1e+5,0,3e+5]);
    set(gca,'XTick',[1e-1,1e+2,1e+5]);
    set(gca,'XTicklabel',{'10^{-1}';'10^{2}';'10^{5}'});
end
set(gca,'YTick',[1e+5,2e+5,3e+5]);
set(gca,'XScale','log');
xlabel('Amino acids (\muM)');
ylabel({'Amino acids';'syn./deg. rates'});


%---------------------------------
%   plot local sensitivity
%---------------------------------

subplot(2,1,2);
hold on;

%   color code
colorPalettes = [255,109,109;255,219,107;72,255,167;61,191,255;198,73,255]/255;

Flag = 0;

deriv               = fnder(spline(ppGppList,AminoAcid(Flag+1,:)),1);
localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./AminoAcid(Flag+1,:));

plot(ppGppList,localSensitivity,'-','Color',colorPalettes(Flag+1,:),'LineWidth',LW);
% plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],...
%     [0,40],'k--','Color',colorPalettes(Flag+1,:));
axis square;
box on;
xlabel('ppGpp (\muM)');
ylabel({'Local sensitivity of';'amino acid pool size'});
axis([40 110 0 40]);
set(gca,'XTick',[50 100]);
set(gca,'YTick',[0,20,40]);


end