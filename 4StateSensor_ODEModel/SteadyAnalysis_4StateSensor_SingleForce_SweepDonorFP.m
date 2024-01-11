close all; clear all; clc;

StudyName = 'SweepDonorFP';
PlotOutDir = StudyName;
mkdir(PlotOutDir);

%% Kinetic Model Params
kT=4.114; % [pN/nm] Approximate value of kT at 298 K

% Mech Switch Rate Constant for Donor
% kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D);
F_half_sweep = [1 3 5 8 25];
m_sweep = [7.6 2.5 1.5 1 1];
SweepPropsFPs.F_half_D = [F_half_sweep F_half_sweep F_half_sweep F_half_sweep F_half_sweep];
SweepPropsFPs.m_D = [m_sweep m_sweep m_sweep m_sweep m_sweep];
SweepPropsFPs.kms0_D = [1 1 1 1 1 .1 .1 .1 .1 .1 .316 .316 .316 .316 .316 3.16 3.16 3.16 3.16 3.16 10 10 10 10 10];
SweepPropsFPs.dxms_D = [.23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23 .23];

% Mech Switch Rate Constant for Acceptor
% kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A);
SweepPropsFPs.F_half_A = [6]*ones(1,25);
SweepPropsFPs.m_A = [1]*ones(1,25);
SweepPropsFPs.kms0_A = [0]*ones(1,25);
SweepPropsFPs.dxms_A = [.23]*ones(1,25);

% MTS Unbinding/Unloading
Toggle_ForceDepUnbinding = 0;
SweepPropsUnbinding.kunload0 = [1];

% Force Sweep
Sweep_F = [0:.1:30];

%% MTS FRET Eff vs Force Model
% Time-Averaged Signal of Single MTS at Specified Force
load('../TheoreticalESHistogramsForMTS\LUT_OriginalTSMod.mat');
MTS_FRETForceFunction = @(F) max(interp1(LUT_OriginalTSMod.Force_app, LUT_OriginalTSMod.FRETEff, F),0);

%% Perform Computation
nSims = length(SweepPropsFPs.F_half_D)*length(SweepPropsUnbinding.kunload0)*length(Sweep_F);
Data = table();
ii=1;
SimulationSetID = 0;
for aa=1:length(SweepPropsFPs.F_half_D)
SimulationSetID = SimulationSetID + 1;
for bb=1:length(SweepPropsUnbinding.kunload0)
for cc=1:length(Sweep_F)

    %% Compute Steady State Fraction for Each Species

    % Get Force
    F = Sweep_F(cc);
    
    % Get Mech Switch Rate Constant for Donor
    F_half_D = SweepPropsFPs.F_half_D(aa);
    m_D = SweepPropsFPs.m_D(aa);
    kms0_D = SweepPropsFPs.kms0_D(aa);
    dxms_D = SweepPropsFPs.dxms_D(aa);
    kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D,kT);

    % Get Mech Switch Rate Constant for Acceptor
    F_half_A = SweepPropsFPs.F_half_A(aa);
    m_A = SweepPropsFPs.m_A(aa);
    kms0_A = SweepPropsFPs.kms0_A(aa);
    dxms_A = SweepPropsFPs.dxms_A(aa);
    kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A,kT);

    % Get Unbinding/Unloading Rate Constant
    if ~Toggle_ForceDepUnbinding
        kunload0 = SweepPropsUnbinding.kunload0(bb);
        kunload = kunload0;
    elseif Toggle_ForceDepUnbinding
        kunload = kunload_TwoPathBellModel(F,k01,xb1,k02,xb2,kT);
    end

    % Compute Steady State Prob of Each Species using Analytical Expressions
    nD1A1 = kunload/(kmsA + kmsD + kunload);
    nD0A1 = (kmsD*kunload)/((kmsA + kunload)*(kmsA + kmsD + kunload));
    nD1A0 = (kmsA*kunload)/((kmsD + kunload)*(kmsA + kmsD + kunload));
    nD0A0 = (kmsA*kmsD*(kmsA + kmsD + 2*kunload))/((kmsA + kunload)*(kmsD + kunload)*(kmsA + kmsD + kunload));

    %% Compute Eapp and Sapp 
    % Use analytical expression that requires all sensors in ensemble have the
    % same force and assumes that D1A1>0

    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 

    % Compute Eapp
    Eapp = E0./(1+nD1A0./nD1A1);

    % Compute Sapp
    S0 = 0.5;
    Sapp = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);

    %% Record Data
    tau0 = 1./kunload;
    Data(ii,:) = table(SimulationSetID,...
        F,kunload0,tau0,...
        F_half_D,m_D,kms0_D,dxms_D,...
        F_half_A,m_A,kms0_A,dxms_A,...
        nD1A1,nD0A0,nD1A0,nD0A0,...
        E0,S0,...
        Eapp,Sapp);
    disp(['Status: ' num2str(round(100*ii/nSims,2)) '% complete.']);
    ii=ii+1;

end
end
end

save(fullfile(PlotOutDir,'AllData_Steady.mat'), 'Data');

%% Figures
load('../TheoreticalESHistogramsForMTS\ESCurves_TSMod_Tension.mat');
set(0,'defaultAxesFontSize',16);

%% [1] Sweep F_half @ kms0=1

% SimulationSetID_toPlot=1:length(unique(Data.SimulationSetID)); % All
SimulationSetID_toPlot = [1:5]; % Sweep F_half @ kms0=1

HighlightPoints_Force = [0 3 6 9 12 15 30];

f = figure;
% Plot ES Curves
for ff=[1 4 7 10 13 16] % F=0,3,6,9,12,15
    p = plot(ESCurves{ff}.Eapp,ESCurves{ff}.Sapp,'k--');
    hold on;
end
% S0 Line
plot([0 .35], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], 'k--');
% Unloaded E0 Line
plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], 'k--');
xlim([0 .3]);
ylim([0 1]);
xticks([0:.1:.3]);
yticks([0:.25:1]);
xlabel('E_{app}');
ylabel('S_{app}');

set(gca,'ColorOrderIndex',1)
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    Eapp=thisData.Eapp;
    Sapp=thisData.Sapp;
    p = plot(Eapp,Sapp);
    thisColor = p.Color;
    
    % Highlight a subset of force points
    for F = HighlightPoints_Force
        Eapp=thisData.Eapp(thisData.F==F);
        Sapp=thisData.Sapp(thisData.F==F);
        % p = plot(Eapp,Sapp,'o','MarkerSize',4);
        % p.Color = thisColor;
        s = scatter(Eapp,Sapp,25,thisColor,'filled','o'); 
    end

end

% Save full figure panel in specified file formats
FigName = 'ES_DonOnly_SweepFhalf_kms0=1';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Plot Donor FP Mech Switch Forward Rate Constants for Each Simulation Set
Toggle_koffref = 1;
f = figure;
F = [0:.1:30];
legend_text = {};
ii=1;
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    
    F_half = thisData(1,:).F_half_D;
    m = thisData(1,:).m_D;
    kms0 = thisData(1,:).kms0_D;
    dxms = thisData(1,:).dxms_D;
    
    if Toggle_koffref
        koff_ref = @(F) kms0.*exp(F.*dxms/kT);
        plot(F,koff_ref(F),'--k'); hold on;
    end
    set(gca,'ColorOrderIndex',ii);
    plot(F,kms(F,F_half,m,kms0,dxms,kT)); hold on;
    
    ii=ii+1;
    
end
xlabel('F [pN]');
ylabel('kms [1/s]'); ylim([0 6]);

% Save full figure panel in specified file formats
FigName = 'RateConstants_DonOnly_SweepFhalf_kms0=1';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% [2] Sweep kms0 @ F_half=6

% SimulationSetID_toPlot=1:length(unique(Data.SimulationSetID)); % All
SimulationSetID_toPlot = [8 13 3 18 23]; % Sweep kms0 @ F_half=6

HighlightPoints_Force = [0 3 6 9 12 15 30];

f = figure;
% Plot ES Curves
for ff=[1 4 7 10 13 16] % F=0,3,6,9,12,15
    p = plot(ESCurves{ff}.Eapp,ESCurves{ff}.Sapp,'k--');
    hold on;
end
% S0 Line
plot([0 .35], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], 'k--');
% Unloaded E0 Line
plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], 'k--');
xlim([0 .3]);
ylim([0 1]);
xticks([0:.1:.3]);
yticks([0:.25:1]);
xlabel('E_{app}');
ylabel('S_{app}');

set(gca,'ColorOrderIndex',1)
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    Eapp=thisData.Eapp;
    Sapp=thisData.Sapp;
    p = plot(Eapp,Sapp);
    thisColor = p.Color;
    
    % Highlight a subset of force points
    for F = HighlightPoints_Force
        Eapp=thisData.Eapp(thisData.F==F);
        Sapp=thisData.Sapp(thisData.F==F);
        % p = plot(Eapp,Sapp,'o','MarkerSize',4);
        % p.Color = thisColor;
        s = scatter(Eapp,Sapp,25,thisColor,'filled','o');
    end

end

% Save full figure panel in specified file formats
FigName = 'ES_DonOnly_Sweepkms0_Fhalf=6';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Plot Donor FP Mech Switch Forward Rate Constants for Each Simulation Set
Toggle_koffref = 1;
f = figure;
F = [0:.1:20];
legend_text = {};
ii=1;
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    
    F_half = thisData(1,:).F_half_D;
    m = thisData(1,:).m_D;
    kms0 = thisData(1,:).kms0_D;
    dxms = thisData(1,:).dxms_D;
    
    if Toggle_koffref
        koff_ref = @(F) kms0.*exp(F.*dxms/kT);
        plot(F,koff_ref(F),'--k'); hold on;
    end
    set(gca,'ColorOrderIndex',ii);
    plot(F,kms(F,F_half,m,kms0,dxms,kT)); hold on;
    
    ii=ii+1;
    
end
xlabel('F [pN]');
ylabel('kms [1/s]'); ylim([0 40]);

% Save full figure panel in specified file formats
FigName = 'RateConstants_DonOnly_Sweepkms0_Fhalf=6';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% Plot Donor FP Mech Switch Forward Rate Constants for Each Simulation Set
SimulationSetID_toPlot=1:length(unique(Data.SimulationSetID)); % All
figure;
F = [0:.1:40];
legend_text = {};
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    
    F_half = thisData(1,:).F_half_D;
    m = thisData(1,:).m_D;
    kms0 = thisData(1,:).kms0_D;
    dxms = thisData(1,:).dxms_D;
    
    plot(F,kms(F,F_half,m,kms0,dxms,kT)); hold on;
    
    legend_text = [legend_text, [num2str(SimulationSetID)]];
    
end
xlabel('F [pN]');
ylabel('kms [1/s]'); ylim([0 10]);
legend(legend_text);

 %% Support Functions
% FP Mechanial Switching, Forward Rate Constant
function out = kms(F,F_half,m,kms0,dxms,kT)
    out = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
end

% Protein Unbinding/Unloading
function out = kunload_TwoPathBellModel(F,k01,xb1,k02,xb2,kT)
    out = k01*exp(xb1*F/kT) + k02*exp(xb2*F/kT);
end
