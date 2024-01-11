clear all;
close all;
clc;

DefineForceDepRateConstantModels;
BondModels{1} = IdealBond_TwoPathBellModelParams;
BondModels{2} = SlipBond_TwoPathBellModelParams;
BondModels{3} = CatchSlipBond_TwoPathBellModelParams;

BondPlotStyles = {'-k','-r','-b'};

Labels = {'Ideal', 'Slip', 'Catch-Slip'};

set(0,'defaultAxesFontSize',8);
set(0,'defaultAxesFontName','Arial');

figure;
dF = .01;
Fmax = 50;
F = [0:dF:Fmax];
for ii=1:length(BondModels)
    
    k01 = BondModels{ii}.k01;
    x1 = BondModels{ii}.xb1;
    k02 = BondModels{ii}.k02;
    x2 = BondModels{ii}.xb2;
    koff{ii} = k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT);
    MeanLifetime{ii} = 1./koff{ii};
    
    plot(F,MeanLifetime{ii},BondPlotStyles{ii});
    hold on;
    
end

xlim([0 35]);
ylim([1E-2 1E1]);
set(gca,'YScale','log');
ylabel('Lifetime [s]');
xlabel('Force [pN]');
legend(Labels);
% grid on;
box on;
title('Force-Dependent Bond Lifetimes');

% Save full figure panel in specified file formats
set(gcf, 'Position', [100   100   400   300]);
Resolution_FigurePanel_png = 300; % Figure Print Resolution for PNG
print(gcf,'BondModelDefinitions.png','-dpng',['-r' num2str(Resolution_FigurePanel_png)],'-painters');
% print(gcf,'BondModelDefinitions_FricValidation.eps','-depsc','-painters');
print(gcf,'BondModelDefinitions.svg','-dsvg','-painters');