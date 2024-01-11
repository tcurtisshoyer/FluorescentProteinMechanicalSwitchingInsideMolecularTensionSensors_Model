close all;
clear all;
clc;

set(0,'defaultAxesFontSize',8);
set(0,'defaultAxesFontName','Arial');

PlotOutDir = 'Plots';
if ~exist(PlotOutDir)
    mkdir(PlotOutDir)
end

Resolution_png = 900; % Figure Print Resolution for PNG

%%  Parameters
kT=4.114; % [pN/nm] Approximate value of kT at 298 K

% FP Mech Switch (Forward Transition)
% koff_combined = P_deltaAlpha(F) * k_MechSwitch(F)
% Base Param Values
deltaG0 = 22*kT; % Dietz and Rief PNAS 2004 (GFP-->GFP_deltaAlpha)
deltaL = 3.2; % [nm] Dietz and Rief PNAS 2004 (Contour Lenght Increase for GFP-->GFP_deltaAlpha)
F_half_base = deltaG0/deltaL;
m_base = 1; % m_base = deltaL/kT;
kms0_base = .33; % [1/s] Ganim and Rief PNAS 2017
kms0 = .33; % [1/s] Ganim and Rief PNAS 2017
dxms = 0.23; % [nm] Ganim and Rief PNAS 2017
dxms_base = 0.23; % [nm] Ganim and Rief PNAS 2017

F_half = F_half_base;
m = m_base;
kms0 = kms0_base;
dxms = dxms_base; 
kun = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);

% Protein Unbinding (Reverse Transition)
% FP recovers immediately after rebinding
kunbind0_base = 1;
kunbind = @(F) kunbind0; % Ideal Bond

% Probability FP Mech Switch
Prob_FPMechSwitch = @(F) kun(F)./(kun(F) + kunbind(F));

%% Parameter Sweep
Data = table();
ii=1;
nSims = 5*5*41*31;
F_half = F_half_base;
m = m_base;
kms0 = kms0_base;
dxms = dxms_base; 
for F_half = logspace(-1,1,5)*F_half_base
for kms0 = flip(logspace(-1,1,5)*kms0_base)
        for F = linspace(0,40,41)
        for kunbind0 = logspace(-2,2,31)*kunbind0_base
            % kun = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
            % kunbind = @(F) kunbind0;
            % Prob_FPMechSwitch = @(F) kun(F)./(kun(F) + kunbind(F));
            kun = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
            kunbind = kunbind0;
            Prob_FPMechSwitch = kun./(kun + kunbind);
            tau = 1/kunbind0;
            Data(ii,:) = table(F,kunbind0,tau,F_half,m,kms0,dxms,kun,kunbind,Prob_FPMechSwitch);
            disp(['Status: ' num2str(round(100*ii/nSims,2)) '% complete.']);
            ii=ii+1;
        end
        end
end
end

save('SingleFP_SweepLoadParams.mat','Data');

%% Make Heatmap Plots of Results
% HeatmapChart Properties: https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.heatmapchart-properties.html
f = figure;
jj=1;
for F_half = logspace(-1,1,5)*F_half_base
for kms0 = logspace(-1,1,5)*kms0_base
    thisData = Data(Data.F_half==F_half & Data.kms0==kms0,:);
    subplot(5,5,jj);
    h = heatmap(thisData,'F','tau','ColorVariable','Prob_FPMechSwitch'); 
    h.Title = ['Fhalf=' num2str(round(F_half,1)) ', kms0=' num2str(round(kms0,2))];
    h.XLabel = 'F [pN]';
    h.YLabel = 'T [s]';
    h.ColorMethod = 'mean';
    h.ColorLimits = [0 1];
    h.Colormap = jet;
    h.GridVisible = 'off';
    h.CellLabelFormat = '%.1f';
    h.NodeChildren(3).YDir='normal';  
    h.FontSize = 4;
    jj=jj+1;
end
end

f.Position = [100 100 1800 950];
FigName = 'ParamSweep_All';
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');
  
%% Make Heatmap Plots of Results
% HeatmapChart Properties: https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.heatmapchart-properties.html
f = figure;
jj=1;
for F_half = logspace(-1,.5,4)*F_half_base
for kms0 = logspace(-1,1,3)*kms0_base
    thisData = Data(Data.F_half==F_half & Data.kms0==kms0,:);
    subplot(4,3,jj);
    h = heatmap(thisData,'F','tau','ColorVariable','Prob_FPMechSwitch'); 
    h.Title = ['Fhalf=' num2str(round(F_half,1)) ', kms0=' num2str(round(kms0,2))];
    h.XLabel = 'F [pN]';
    h.YLabel = 'T [s]';
    h.ColorMethod = 'mean';
    h.ColorLimits = [0 1];
    h.Colormap = jet;
    % h.ColorbarVisible = 'off';
    h.GridVisible = 'off';
    h.CellLabelFormat = '%.1f';
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    h.NodeChildren(3).YDir='normal'; 
    h.FontSize = 8;
    jj=jj+1;
end
end

f.Position = [100 100 1200 950];
FigName = 'ParamSweep_Subset';
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% Schematic Plots: FP Params - Generic

% Param Sweeps
F = [0:.1:40];

% Sweep F_half around base value
F_half = 20;
m = m_base;
kms0 = 1;
dxms = dxms_base; 
f = figure;
legends = {};
for F_half = [10 20 30]    
    koff_combined = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
    plot(F,koff_combined(F)); hold on;
    % legends = [legends; ['Fhalf = ' num2str(round(F_half,1))]];
end 
koff_ref = @(F) kms0.*exp(F.*dxms/kT);
plot(F,koff_ref(F),'--k');
% legend(legends);
ylabel('k [1/s]');
xlabel('F [pN]');
% set(gca, 'YScale', 'log'); ylim([.01 100]);
ylim([0 10]); 
xlim([0 40]);
% Save figure
FigName = 'ParamDemo_Sweep_Fhalf';
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Sweep kms0 around base value
F_half = 20;
m = m_base;
kms0 = 1;
dxms = dxms_base; 
f = figure;
legends = {};
for kms0 = logspace(-.5,.5,3)*1    
    koff_combined = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
    plot(F,koff_combined(F)); hold on;
    % legends = [legends; ['kms0 = ' num2str(round(kms0,1))]];
    koff_ref = @(F) kms0.*exp(F.*dxms/kT);
    plot(F,koff_ref(F),'--k');
end
% legend(legends);
ylabel('k [1/s]');
xlabel('F [pN]');
% set(gca, 'YScale', 'log'); 
ylim([0 30]); 
xlim([0 40]);
% Save figure
FigName = 'ParamDemo_Sweep_kms0';
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');
 
%% Schematic Plots: FP Params - Actual
if 0
% Param Sweeps
F = [0:.1:40];

% Sweep F_half around base value
F_half = F_half_base;
m = m_base;
kms0 = kms0_base;
dxms = dxms_base; 
figure;
legends = {};
for F_half = logspace(-1,1,5)*F_half_base    
    koff_combined = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
    plot(F,koff_combined(F)); hold on;
    legends = [legends; ['Fhalf = ' num2str(round(F_half,1))]];
end
legend(legends);
ylabel('k [1/s]');
xlabel('F [pN]');
% set(gca, 'YScale', 'log'); ylim([.01 100]);
ylim([0 40]); 
xlim([0 40]);

% Sweep kms0 around base value
F_half = F_half_base;
m = m_base;
kms0 = kms0_base;
dxms = dxms_base; 
figure;
legends = {};
for kms0 = logspace(-1,1,5)*kms0_base    
    koff_combined = @(F) (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
    plot(F,koff_combined(F)); hold on;
    legends = [legends; ['kms0 = ' num2str(round(kms0,1))]];
end
legend(legends);
ylabel('k [1/s]');
xlabel('F [pN]');
% set(gca, 'YScale', 'log'); 
ylim([0 40]); 
xlim([0 40]);
end      
        
 %% Support Functions
% FP Mechanial Switching, Forward Rate Constant
function out = kms(F,F_half,m,kms0,dxms,kT)
    out = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
end

% Protein Unbinding/Unloading
function out = kunload_TwoPathBellModel(F,k01,xb1,k02,xb2,kT)
    out = k01*exp(xb1*F/kT) + k02*exp(xb2*F/kT);
end