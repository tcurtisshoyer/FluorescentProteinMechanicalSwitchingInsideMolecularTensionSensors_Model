close all; clear all; clc;

StudyName = 'SweepCatchBond';
PlotOutDir = StudyName;
mkdir(PlotOutDir);

% What? 9 Diff Catch-Slip Bonds with Different Fb values

%% Kinetic Model Params
kT=4.114; % [pN/nm] Approximate value of kT at 298 K

% Mech Switch Rate Constant for Donor
% kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D);
SweepPropsFPs.F_half_D = [5 5 5];
SweepPropsFPs.m_D = [1.5 1.5 1.5];
SweepPropsFPs.kms0_D = [0 0 1];
SweepPropsFPs.dxms_D = [.23 .23 .23];

% Mech Switch Rate Constant for Acceptor
% kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A);
SweepPropsFPs.F_half_A = [5 5 5];
SweepPropsFPs.m_A = [1.5 1.5 1.5];
SweepPropsFPs.kms0_A = [0 1 0];
SweepPropsFPs.dxms_A = [.23 .23 .23];

% MTS Unbinding/Unloading
% Specified as Parameters for 2 Pathway Catch-Slip Bond Model: 
% koff(F)k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT)
% Make Catch-Slip Bond Models with Different Fb
id_BondModel = 1;
koff0 = 1; % Intrinsic Off-Rate
for Fb = [Inf 1 3 5 8 25] % Coefficent [pN] in Bell Model
    xb = kT/Fb; % Coefficent [nm] in Bell Model
    CatchSlipBond_TwoPathBellModelParams.k01 = koff0*.9;
    CatchSlipBond_TwoPathBellModelParams.xb1 = -1*2*xb;
    CatchSlipBond_TwoPathBellModelParams.k02 = koff0*.1;
    CatchSlipBond_TwoPathBellModelParams.xb2 = xb;
    BondModels{id_BondModel} = CatchSlipBond_TwoPathBellModelParams;
    id_BondModel = id_BondModel+1;
end
SweepPropsUnbinding.bondModelID = [1 2 3 4 5 6]; % Specifies which bond model to use
SweepPropsUnbinding.intrinsicRateConstantScaling = [1 1 1 1 1 1]; % Specifies scaling of koff

% Force Sweep
Sweep_F = [0:.1:30];

%% MTS FRET Eff vs Force Model
% Time-Averaged Signal of Single MTS at Specified Force
load('../../TheoreticalESHistogramsForMTS\LUT_OriginalTSMod.mat');
MTS_FRETForceFunction = @(F) max(interp1(LUT_OriginalTSMod.Force_app, LUT_OriginalTSMod.FRETEff, F),0);

%% Perform Computation
nSims = length(SweepPropsFPs.F_half_D)*length(SweepPropsUnbinding.bondModelID)*length(Sweep_F);
Data = table();
ii=1;
SimulationSetID = 0;
for aa=1:length(SweepPropsFPs.F_half_D)
for bb=1:length(SweepPropsUnbinding.bondModelID)
    SimulationSetID = SimulationSetID + 1;
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
    thisBondModel = BondModels{SweepPropsUnbinding.bondModelID(bb)};
    thisIntrinsicRateConstantScaling = SweepPropsUnbinding.intrinsicRateConstantScaling(bb);
    k01 = thisBondModel.k01 * thisIntrinsicRateConstantScaling;
    xb1 = thisBondModel.xb1;
    k02 = thisBondModel.k02 * thisIntrinsicRateConstantScaling;
    xb2 = thisBondModel.xb2;
    kunload = k01.*exp(xb1.*F/kT) + k02.*exp(xb2.*F/kT);

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
    tau = 1./kunload;
    if SweepPropsUnbinding.bondModelID(bb)==1
        BondType = {'Ideal'};
    else
        BondType = {'CatchSlip'};
    end
    kunload_k01 = thisBondModel.k01;
    kunload_xb1 = thisBondModel.xb1;
    kunload_k02 = thisBondModel.k02;
    kunload_xb2 = thisBondModel.xb2;
    Data(ii,:) = table(SimulationSetID,...
        F,...
        BondType,...
        kunload,tau,...
        kunload_k01,kunload_xb1,kunload_k02,kunload_xb2,...
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
load(fullfile(PlotOutDir,'AllData_Steady.mat'));
load('../../TheoreticalESHistogramsForMTS\ESCurves_TSMod_Tension.mat');
set(0,'defaultAxesFontSize',16);

%% [1] Acceptor Only: Ideal, 5x Catch-Slip Bonds with Different Fb values

% SimulationSetID_toPlot=1:length(unique(Data.SimulationSetID)); % All
SimulationSetID_toPlot = [7:12];

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
FigName = 'ES_Acceptor';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Plot Unbinding Rate Constant for Each Simulation Set
f = figure;
F = [0:.1:30];
legend_text = {};
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    
    k01 = thisData(1,:).kunload_k01;
    xb1 = thisData(1,:).kunload_xb1;
    k02 = thisData(1,:).kunload_k02;
    xb2 = thisData(1,:).kunload_xb2;
    
    kunload = k01.*exp(xb1.*F/kT) + k02.*exp(xb2.*F/kT);
    
    plot(F,kunload); hold on;
    
    legend_text = [legend_text, ['Fb = ' num2str(round(kT/xb2,2))]];
    
end

% Add plot of FP Mech Switch Param Overlaid (assume same for all sims)
F_half = thisData(1,:).F_half_A;
m = thisData(1,:).m_A;
kms0 = thisData(1,:).kms0_A;
dxms = thisData(1,:).dxms_A;
plot(F,kms(F,F_half,m,kms0,dxms,kT),'--k'); hold on;
legend_text = [legend_text, 'FP Mech Switch'];

xlabel('F [pN]');
ylabel('k [1/s]'); ylim([0 10]);
legend(legend_text);

% Save full figure panel in specified file formats
FigName = 'UnbindingRateConstants_Acceptor_LinLin';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Save full figure panel in specified file formats
set(gca, 'YScale', 'log');
ylim([10^(-2) 10^(2)]);
FigName = 'UnbindingRateConstants_Acceptor_LinLog';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% [2] Donor Only: Ideal, 5x Catch-Slip Bonds with Different Fb values

% SimulationSetID_toPlot=1:length(unique(Data.SimulationSetID)); % All
SimulationSetID_toPlot = [13:18];

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
FigName = 'ES_Donor';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

% Plot Unbinding Rate Constant for Each Simulation Set
f = figure;
F = [0:.1:30];
legend_text = {};
for SimulationSetID=SimulationSetID_toPlot
    thisData = Data(Data.SimulationSetID==SimulationSetID,:);
    
    k01 = thisData(1,:).kunload_k01;
    xb1 = thisData(1,:).kunload_xb1;
    k02 = thisData(1,:).kunload_k02;
    xb2 = thisData(1,:).kunload_xb2;
    
    kunload = k01.*exp(xb1.*F/kT) + k02.*exp(xb2.*F/kT);
    
    plot(F,kunload); hold on;
    
    legend_text = [legend_text, ['Fb = ' num2str(round(kT/xb2,2))]];
    
end

% Add plot of FP Mech Switch Param Overlaid (assume same for all sims)
F_half = thisData(1,:).F_half_D;
m = thisData(1,:).m_D;
kms0 = thisData(1,:).kms0_D;
dxms = thisData(1,:).dxms_D;
plot(F,kms(F,F_half,m,kms0,dxms,kT),'--k'); hold on;
legend_text = [legend_text, 'FP Mech Switch'];

xlabel('F [pN]');
ylabel('k [1/s]'); ylim([0 10]);
legend(legend_text);

% Save full figure panel in specified file formats
FigName = 'UnbindingRateConstants_Donor';
Resolution_png = 300;
print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

%% Support Functions
% FP Mechanial Switching, Forward Rate Constant
function out = kms(F,F_half,m,kms0,dxms,kT)
    out = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
end

% Protein Unbinding/Unloading
function out = kunload_TwoPathBellModel(F,k01,xb1,k02,xb2,kT)
    out = k01*exp(xb1*F/kT) + k02*exp(xb2*F/kT);
end
