close all; clear all; clc;

% ES Curves for TSMod at Different Forces

% ES Curve = Load Isocline: curve of constant molecular load (and thus E0)
% with varying amounts of broken A only or broken D only.

% MTS FRET Eff vs Force Model
% Time-Averaged Signal of Single MTS at Specified Force
load('LUT_OriginalTSMod.mat');
MTS_FRETForceFunction = @(F) max(interp1(LUT_OriginalTSMod.Force_app, LUT_OriginalTSMod.FRETEff, F),0);

nSensor = 1;

ForceSweep = [0:1:15];
for ff=1:length(ForceSweep)
    thisESCurve = table();
    F = ForceSweep(ff);
    
    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 
    
    S0 = 0.5;

    % Only Acceptor Breaks (--> Free Donor)
    fractionDBroken = 0;
    for fractionABroken = [0:.001:1]
        nD0A1 = nSensor*fractionDBroken; % Only Donor Broken
        nD1A0 = nSensor*fractionABroken; % Only Acceptor Broken
        nD0A0 = 0; % Both Donor and Acceptor Broken
        nD1A1 = nSensor-nD0A1-nD1A0-nD0A0; % Ideal Sensor 
        if nD1A1==0 % Eapp is undefined for 0 Ideal Sensors. Set to Eapp for very small # Ideal Sensors.
            nD1A1=1E-5;
        end
        Sapp = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1); % Apparent Stoichiometry
        Eapp = E0./(1+nD1A0./nD1A1); % Apparent FRET Eff
        thisESCurve = [thisESCurve; table(F, E0, S0, Eapp, Sapp, fractionDBroken, fractionABroken,...
            nSensor, nD1A1, nD1A0, nD0A1, nD0A0)];
    end

    % Only Donor Breaks (--> Free Acceptor)
    fractionABroken = 0;
    for fractionDBroken = [.001:.001:1]
        nD0A1 = nSensor*fractionDBroken; % Only Donor Broken
        nD1A0 = nSensor*fractionABroken; % Only Acceptor Broken
        nD0A0 = 0; % Both Donor and Acceptor Broken
        nD1A1 = nSensor-nD0A1-nD1A0-nD0A0; % Ideal Sensor 
        if nD1A1==0 % Eapp is undefined for 0 Ideal Sensors. Set to Eapp for very small # Ideal Sensors.
            nD1A1=1E-5;
        end
        Sapp = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1); % Apparent Stoichiometry
        Eapp = E0./(1+nD1A0./nD1A1); % Apparent FRET Eff
        thisESCurve = [thisESCurve; table(F, E0, S0, Eapp, Sapp, fractionDBroken, fractionABroken,...
            nSensor, nD1A1, nD1A0, nD0A1, nD0A0)];
   end
    
    % Sort by Ascending Sapp
    thisESCurve = sortrows(thisESCurve,'Sapp');
    
    ESCurves{ff} = thisESCurve;
end

ESCurve_0pN = ESCurves{1};
ESCurve_1pN = ESCurves{2};
ESCurve_2pN = ESCurves{3};
ESCurve_3pN = ESCurves{4};
ESCurve_4pN = ESCurves{5};
ESCurve_5pN = ESCurves{6};
ESCurve_6pN = ESCurves{7};

save('ESCurves_TSMod_Tension.mat',...
    'ESCurves',...
    'ESCurve_0pN',...
    'ESCurve_1pN',...
    'ESCurve_2pN',...
    'ESCurve_3pN',...
    'ESCurve_4pN',...
    'ESCurve_5pN',...
    'ESCurve_6pN');

%% Theoretical ES Curves for F=0,3,6,9,12,15
set(0,'defaultAxesFontSize',10);
set(0,'defaultAxesFontName','Arial');

f=figure;
leg={};
ii=1;
for ff=[1 4 7 10 13 16] % F=0,3,6,9,12,15
    p(ii) = plot(ESCurves{ff}.Eapp,ESCurves{ff}.Sapp);
    hold on;
    leg=[leg, ['F = ' num2str(ForceSweep(ff)) ' pN'] ];
    ii=ii+1;
end
% S0 Line
plot([0 .35], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], 'k--');
% Unloaded E0 Line
plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], 'k--');
xlim([0 .35]);
ylim([0 1]);
xticks([0:.1:.3]);
yticks([0:.25:1]);
xlabel('E_{app}');
ylabel('S_{app}');
legend(leg);

% Save full figure panel in specified file formats
PlotOutDir = 'Plots';
FigName = 'TheoreticalESCurves';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% E0 Points for F=0,3,6,9,12,15
f = figure;
leg={};
ii=1;
for ff=[1 4 7 10 13 16] % F=0,3,6,9,12,15
    fractionABroken = 0;
    fractionDBroken = 0;
    jj = find(thisESCurve.fractionABroken==fractionABroken & thisESCurve.fractionDBroken==fractionDBroken);
    s(ii) = scatter(ESCurves{ff}.Eapp(jj),ESCurves{ff}.Sapp(jj));
    s(ii).Marker = 'o';
    s(ii).SizeData = 50;
    s(ii).MarkerFaceColor = get(p(ii),'Color');
    s(ii).MarkerEdgeColor = get(p(ii),'Color');
    hold on;
    leg=[leg, ['F = ' num2str(ForceSweep(ff)) ' pN'] ];
    ii=ii+1;
end
% S0 Line
plot([0 .35], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], 'k--');
% Unloaded E0 Line
plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], 'k--');
xlim([0 .35]);
ylim([0 1]);
xticks([0:.1:.3]);
yticks([0:.25:1]);
xlabel('E_{app}');
ylabel('S_{app}');
legend(leg);
box on;

FigName = 'E0PointsOnly';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% One ES Curve with Inicators of Fraction Broken
plotForce = 3;
thisESCurve = ESCurves{find(ForceSweep==plotForce)};
f = figure;
plot(thisESCurve.Eapp,thisESCurve.Sapp,'k');
hold on;
xlabel('E_{app}');
ylabel('S_{app}');
xlim([0 .35]);
ylim([0 1]);
xticks([0:.1:.3]);
yticks([0:.25:1]);
% S0 Line
plot([0 .35], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], 'k--');
% Unloaded E0 Line
plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], 'k--');

% Plot some points for broken acceptor
fractionDBroken = 0;
for fractionABroken = [.1 .5 .9]
    jj = find(thisESCurve.fractionABroken==fractionABroken & thisESCurve.fractionDBroken==fractionDBroken);
    s = scatter(thisESCurve.Eapp(jj),thisESCurve.Sapp(jj));
    s.Marker = 'o';
    s.SizeData = 20+200*fractionABroken;
    s.MarkerFaceColor = 'y';
    s.MarkerEdgeColor = 'k';
end

% Plot some points for broken donor
fractionABroken = 0;
for fractionDBroken = [.1 .5 .9]
    jj = find(thisESCurve.fractionABroken==fractionABroken & thisESCurve.fractionDBroken==fractionDBroken);
    s = scatter(thisESCurve.Eapp(jj),thisESCurve.Sapp(jj));
    s.Marker = 'o';
    s.SizeData = 20+200*fractionDBroken;
    s.MarkerFaceColor = 'c';
    s.MarkerEdgeColor = 'k';
end

FigName = 'ESCurve_3pN_wFractionBroken';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% Fret Eff vs Force Curve with Points for F=0,3,6,9,12,15
f = figure;
F = [0:.01:30];
plot(F,MTS_FRETForceFunction(F),'k'); hold on;
hold on;
ylabel('E_0');
xlabel('Force [pN]');
for ff=[1 4 7 10 13 16] % F=0,3,6,9,12,15
    s = scatter(ForceSweep(ff),MTS_FRETForceFunction(ForceSweep(ff)));
    s.Marker = 'o';
    s.SizeData = 50;
    s.MarkerFaceColor = get(p(ff),'Color');
    s.MarkerEdgeColor = get(p(ff),'Color');
    hold on;
end

FigName = 'EffForceCalibrationCurve_TSMod';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');