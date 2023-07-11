function plot_sensitivity_values_vary_ppGpp(ppGppList,AminoAcid,Ribosome,growthRate,...
    ppGppSynRate,ppGppDegRate,IndexLowerLimit,IndexUpperLimit,IndexMaxGrowth,hostparams, Flag)

%   This parameter is used to remove boundary effects of 
offset = 5;

figure();
hold on;

LW = 2;

%   ribosome sensitivity
deriv               = fnder(spline(ppGppList,Ribosome(Flag+1,:)),1);
localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./Ribosome(Flag+1,:));

subplot(3,1,1);
hold on;
plot(ppGppList(IndexLowerLimit+offset:end),localSensitivity(IndexLowerLimit+offset:end),'k-','LineWidth',LW);
plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,10],'k--');
axis square;
box on;
xlabel('ppGpp (\muM)');
ylabel({'Sensitivity of';'ribosome concentration'});
axis([0 200 0 10]);
set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0,5,10]);

%   amino acid sensitivity
deriv               = fnder(spline(ppGppList,AminoAcid(Flag+1,:)),1);
localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./AminoAcid(Flag+1,:));

subplot(3,1,2);
hold on;
plot(ppGppList(IndexLowerLimit+offset:end),localSensitivity(IndexLowerLimit+offset:end),'k-','LineWidth',LW);
plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,40],'k--');
axis square;
box on;
xlabel('ppGpp (\muM)');
ylabel({'Sensitivity of';'amino acid concentration'});
axis([0 200 0 40]);
set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0,20,40]);

%   sensitivity of ppGpp synthesis rate
deriv               = fnder(spline(ppGppList,ppGppSynRate(Flag+1,:)),1);
localSensitivity    = abs(ppval(deriv,ppGppList).*ppGppList./ppGppSynRate(Flag+1,:));

subplot(3,1,3);
hold on;
plot(ppGppList(IndexLowerLimit+offset:end),localSensitivity(IndexLowerLimit+offset:end),'k-','LineWidth',LW);
plot([ppGppList(IndexMaxGrowth(Flag+1)),ppGppList(IndexMaxGrowth(Flag+1))],[0,40],'k--');
axis square;
box on;
xlabel('ppGpp (\muM)');
ylabel({'Sensitivity of';'ribosome concentration'});
axis([0 200 0 40]);
set(gca,'XTick',[0 100 200]);
set(gca,'YTick',[0,20,40]);


end

