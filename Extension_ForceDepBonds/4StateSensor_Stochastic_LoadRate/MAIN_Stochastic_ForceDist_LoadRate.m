close all; clear all; clc;

% FP Mech Switch Params are 4 Base Cases: None, A only, D only, Both
% Sweep Loading Rate
% MTS Unbinding is Ideal, Slip, or Catch-Slip Bond (2 Path Bell Models)
% Imaging Noise Cases: Off

StudyName = 'LoadingRateDependence_01';
PlotOutDir = StudyName;
mkdir(PlotOutDir);
AllPixelDataOutDir = fullfile('C:\Users\tcsho\Desktop\SimData_FPMechSwitchModel_ForceDepBonds',StudyName);
mkdir(AllPixelDataOutDir);

Toggle_OutputTimeCourseOfEnsemble = 1;
Toggle_OutputLogOfAllTransitions = 1;

%% Simulation Algorithm Parameters
SimulationParams.tend = 100; % Total simulation time [s]
SimulationParams.interval_LogEnsembleData = .01; % Time interval for logging sensor data 
SimulationParams.timeSpan_LogTransitionData = [0.5 1]*SimulationParams.tend; % Time span over which to log linkage transisitions
SimulationParams.maxAllowedTimeStep = 1; % Max Allowed Timestep in Gillespie SSA. Meaningful for dynamic loading profiles ( i.e. F=f(t) ).If dt larger than this value, advance time and update mechanics without executing a sensor transition

%% Model Size Parameters
nSensors = 50; % Used this for previous Same Force sims: 100
nPixelsToSimulate = 100; % Used this for previous Same Force sims: 1000

%% Kinetic Model Params
kT=4.114; % [pN/nm] Approximate value of kT at 298 K

% Mech Switch Rate Constant for Donor
% kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D);
SweepPropsFPs.F_half_D = [5 5 5 5];
SweepPropsFPs.m_D = [1.5 1.5 1.5 1.5];
SweepPropsFPs.kms0_D = [0 0 1 1];
SweepPropsFPs.dxms_D = [.23 .23 .23 .23];

% Mech Switch Rate Constant for Acceptor
% kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A);
SweepPropsFPs.F_half_A = [5 5 5 5];
SweepPropsFPs.m_A = [1.5 1.5 1.5 1.5];
SweepPropsFPs.kms0_A = [0 1 0 1];
SweepPropsFPs.dxms_A = [.23 .23 .23 .23];

% MTS Unbinding/Unloading
% All sensors have same parameters in a given sim ("single value")
% Specified as Parameters for 2 Pathway Catch-Slip Bond Model: 
% koff(F)k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT)
DefineForceDepRateConstantModels; % *See this file for definition of bond models
BondModels{1} = IdealBond_TwoPathBellModelParams;
BondModels{2} = SlipBond_TwoPathBellModelParams;
BondModels{3} = CatchSlipBond_TwoPathBellModelParams;
SweepPropsLoading.bondModelID = [1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3]; % Specifies which bond model to use
SweepPropsLoading.intrinsicRateConstantScaling = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Specifies scaling of koff

% Load Rate Sweep [pN/s]
SweepPropsLoading.LoadingRate = [logspace(-1,2,7),logspace(-1,2,7),logspace(-1,2,7)];

%% MTS FRET Eff vs Force Model
% Time-Averaged Signal of Single MTS at Specified Force
load('../../TheoreticalESHistogramsForMTS\LUT_OriginalTSMod.mat');
MTS_FRETForceFunction = @(F) max(interp1(LUT_OriginalTSMod.Force_app, LUT_OriginalTSMod.FRETEff, F),0);

%% Perform Simulations
nSims = length(SweepPropsFPs.F_half_D)*length(SweepPropsLoading.LoadingRate)*1;
Data_OneRowPerPixel = table();
Data_OneRowPerSimID = table();
SimulationSetID = 0;
SimID=0; 
ii = 1;
for Toggle_Noise = [0] % Simulate Noise in Observed Signal: 0=No Noise, 1=Noise
for aa=1:length(SweepPropsFPs.F_half_D)
for bb=1:length(SweepPropsLoading.LoadingRate)
    
    SimulationSetID = SimulationSetID + 1;
    
    % Get Mech Switch Rate Constant for Donor
    % kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D,kT);
    F_half_D = SweepPropsFPs.F_half_D(aa);
    m_D = SweepPropsFPs.m_D(aa);
    kms0_D = SweepPropsFPs.kms0_D(aa);
    dxms_D = SweepPropsFPs.dxms_D(aa);

    % Get Mech Switch Rate Constant for Acceptor
    % kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A,kT)
    F_half_A = SweepPropsFPs.F_half_A(aa);
    m_A = SweepPropsFPs.m_A(aa);
    kms0_A = SweepPropsFPs.kms0_A(aa);
    dxms_A = SweepPropsFPs.dxms_A(aa);
    
    % Get unbinding rate constant
    % Catch-Slip Bond Model: koff(F)k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT)
    thisBondModelID = SweepPropsLoading.bondModelID(bb);
    thisintrinsicRateConstantScaling = SweepPropsLoading.intrinsicRateConstantScaling(bb);
    SensorUnbindingRateConstant_k01 = BondModels{thisBondModelID}.k01 * thisintrinsicRateConstantScaling;
    SensorUnbindingRateConstant_xb1 = BondModels{thisBondModelID}.xb1;
    SensorUnbindingRateConstant_k02 = BondModels{thisBondModelID}.k02 * thisintrinsicRateConstantScaling;
    SensorUnbindingRateConstant_xb2 = BondModels{thisBondModelID}.xb2;

    % Set loading rate
    LoadingRate = SweepPropsLoading.LoadingRate(bb);
           
    SimID = SimID + 1;
    UniqueID = ['SimSet_' num2str(SimulationSetID) '_SimID_' num2str(SimID)];
    
    %% Simulate Each Pixel
    SSAData = cell(nPixelsToSimulate,1);
    CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState = []; % Stores the Force, Unbinding Rate Constant, Donor State, and Acceptor state for each sensor for all pixels in a given parameter combination 
    for idPixel=1:nPixelsToSimulate
        
        %% Simulate MTS Ensemble w/ Kinetic Model for FP Mech Switching in MTS
        % Number of Sensors
        ModelParams.numSensors = nSensors;
        % Sensor Loading Rate - All sensors have same parameters in a given sim ("single value")
        % Set above in param combo loop
        ModelParams.LoadingRate = LoadingRate;
        % Sensor Unbinding Rate Constant - All sensors have same parameters in a given sim ("single value")
        % Specified as Parameters for 2 Pathway Catch-Slip Bond Model: 
        % koff(F)k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT)
        % Set above in param combo loop
        ModelParams.SensorUnbindingRateConstant_k01 = SensorUnbindingRateConstant_k01;
        ModelParams.SensorUnbindingRateConstant_xb1 = SensorUnbindingRateConstant_xb1;
        ModelParams.SensorUnbindingRateConstant_k02 = SensorUnbindingRateConstant_k02;
        ModelParams.SensorUnbindingRateConstant_xb2 = SensorUnbindingRateConstant_xb2;
        % Mech Switching Paramters for Donor - All sensors have same kinetics
        % kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D);
        ModelParams.F_half_D = F_half_D;
        ModelParams.m_D = m_D;
        ModelParams.kms0_D = kms0_D;
        ModelParams.dxms_D = dxms_D;
        % Mech Switching Paramters for Acceptor - All sensors have same kinetics
        % kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A);
        ModelParams.F_half_A = F_half_A;
        ModelParams.m_A = m_A;
        ModelParams.kms0_A = kms0_A;
        ModelParams.dxms_A = dxms_A;

        % Simulation MTS Ensemble
        [AllSensors_Force,AllSensors_Lifetime,AllSensors_State_Donor,AllSensors_State_Acceptor,num_D1A1,num_D1A0,num_D0A1,num_D0A0,SimulationPerformanceStats,DataLog_Ensemble,DataLog_Transitions,DataLog_LifetimesAtTurnoverOnly,DataLog_AllVariableNames] = ...
            SSA_4SS_LoadRate_ForceDepBonds(ModelParams,SimulationParams,Toggle_OutputTimeCourseOfEnsemble,Toggle_OutputLogOfAllTransitions);
        
        % Add SSA Data Logs to SSAData (to be exported later as single file
        % containing entries for each Pixel ID)
        SSAData{idPixel}.SimulationSetID = SimulationSetID;
        SSAData{idPixel}.SubSimID = SimID;
        SSAData{idPixel}.UniqueID = UniqueID;
        SSAData{idPixel}.idPixel=idPixel;
        SSAData{idPixel}.ModelParams=ModelParams;
        SSAData{idPixel}.SimulationParams=SimulationParams;
        SSAData{idPixel}.SimulationPerformanceStats=SimulationPerformanceStats;
        FinalSensorPopulation_Force_Lifetime_DState_AState = [AllSensors_Force,AllSensors_Lifetime,AllSensors_State_Donor,AllSensors_State_Acceptor];            
        SSAData{idPixel}.FinalSensorPopulation_Force_UnbindRate_DState_AState=FinalSensorPopulation_Force_Lifetime_DState_AState;
        if Toggle_OutputTimeCourseOfEnsemble
            SSAData{idPixel}.DataLog_AllVariableNames=DataLog_AllVariableNames;
            SSAData{idPixel}.DataLog_Ensemble=DataLog_Ensemble;
        end
        if Toggle_OutputLogOfAllTransitions
            SSAData{idPixel}.DataLog_AllVariableNames=DataLog_AllVariableNames;
            SSAData{idPixel}.DataLog_Transitions=DataLog_Transitions;  
        end

        % Add sensor population to CombinedSensorPopulation_idPx_F_DState_AState
        CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState = [CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState;
            repmat(idPixel,[size(FinalSensorPopulation_Force_Lifetime_DState_AState,1) 1]),...
            FinalSensorPopulation_Force_Lifetime_DState_AState];

        % Plot Timetrace of Species Numbers for Five Representative Sims
        if Toggle_OutputTimeCourseOfEnsemble && idPixel <= 5
            f=figure('Visible','off');
            plot(DataLog_Ensemble(:,1),DataLog_Ensemble(:,2:5));
            xlabel('t [s]'); xlim([0 1.1*SimulationParams.tend]);
            ylabel('Species #');
            legend({'D1A1','D1A0','D0A1','D0A0'},'Location','northeastoutside');
            Resolution_Figure_png = 300;
            print(f,fullfile(AllPixelDataOutDir,['SSA_Timecourse_' UniqueID '_pxID_' num2str(idPixel) '.png']),'-dpng',['-r' num2str(Resolution_Figure_png)],'-painters');
            print(f,fullfile(AllPixelDataOutDir,['SSA_Timecourse_' UniqueID '_pxID_' num2str(idPixel) '.svg']),'-dsvg','-painters');
            close(f);
            % For Diagnostics - Mean Sensor Turnover Lifetime in SSA
            % DataLog_TurnoverTransitions = DataLog_Transitions(DataLog_Transitions(:,3)==0,:); % idStateVariable==0 --> Turnover Transition
            % mean_SensorTurnoverLifetime = mean(DataLog_TurnoverTransitions(:,9));
        end
        
        %% Simulate Three Channel FRET Measurements of MTS Ensemble
            
        % Obtain FRET Eff for Ideal Sensor State for Each Sensor based
        % on Sensor Force Distribution
        AllSensors_E0 = MTS_FRETForceFunction(AllSensors_Force);
        AllSensors_S0 = 0.5*ones(size(AllSensors_E0));

        % Store D#A# State for Each Sensor to pass to
        % simulateThreeChannelFRET_E0dist.m
        AllSensors_D1A1 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1; % Set to 1 if D1A1, else 0
        AllSensors_D1A0 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0; % Set to 1 if D1A0, else 0
        AllSensors_D0A1 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1; % Set to 1 if D0A1, else 0
        AllSensors_D0A0 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0; % Set to 1 if D0A0, else 0

        PixelData = simulateThreeChannelFRET_E0dist(AllSensors_D1A1, AllSensors_D1A0, AllSensors_D0A1, AllSensors_D0A0, AllSensors_E0, Toggle_Noise);
        Iaa = PixelData.Iaa; 
        Idd = PixelData.Idd;
        Ida = PixelData.Ida;
        Fc = PixelData.Fc;
        Eapp = PixelData.Eapp;
        Sapp = PixelData.Sapp;
        
        %% Log Data
        Force_mean = mean(AllSensors_Force);
        Force_median = median(AllSensors_Force);
        Force_std = std(AllSensors_Force);
        Force_min = min(AllSensors_Force);
        Force_max = max(AllSensors_Force);

        % Note: This is the current lifetime counter for each sensor. Thus,
        % it is not useful for computing statistics on bound lifetime,
        % which must be computed from the transition log by analyzing the
        % lifetime associated with unbinding transitions.
        Lifetime_mean = mean(AllSensors_Lifetime);
        Lifetime_median = median(AllSensors_Lifetime);
        Lifetime_std = std(AllSensors_Lifetime);
        Lifetime_min = min(AllSensors_Lifetime);
        Lifetime_max = max(AllSensors_Lifetime);
        
        E0_mean = mean(AllSensors_E0);
        E0_median = median(AllSensors_E0);
        E0_std = std(AllSensors_E0);

        S0_mean = mean(AllSensors_S0);
        S0_median = median(AllSensors_S0);
        S0_std = std(AllSensors_S0);

        nD = num_D1A1+num_D1A0;
        nA = num_D1A1+num_D0A1;
        fractionDonorBroken = (nSensors-nD)/nSensors;
        fractionAcceptorBroken = (nSensors-nA)/nSensors;
        fractionSensorsBroken = (nSensors-num_D1A1)/nSensors;

        if SweepPropsLoading.bondModelID(bb)==1
            BondType = {'Ideal'};
        elseif SweepPropsLoading.bondModelID(bb)==2
            BondType = {'Slip'};
        elseif SweepPropsLoading.bondModelID(bb)==3
            BondType = {'CatchSlip'};
        end

        Data_OneRowPerPixel(ii,:) = table(...
            SimulationSetID,SimID,idPixel,...
            Toggle_Noise,...
            LoadingRate,...
            BondType,...
            SensorUnbindingRateConstant_k01,SensorUnbindingRateConstant_xb1,...
            SensorUnbindingRateConstant_k02,SensorUnbindingRateConstant_xb2,...
            F_half_D,m_D,kms0_D,dxms_D,...
            F_half_A,m_A,kms0_A,dxms_A,...
            nSensors,...
            Force_mean,Force_median,Force_std,Force_min,Force_max,...
            Lifetime_mean,Lifetime_median,Lifetime_std,Lifetime_min,Lifetime_max,...
            E0_mean,E0_median,E0_std,...
            S0_mean,S0_median,S0_std,...
            Eapp,Sapp,...
            num_D1A1,num_D0A1,num_D1A0,num_D0A0,...
            fractionDonorBroken,fractionAcceptorBroken,fractionSensorsBroken,Iaa,Idd,Ida,Fc);
        
        ii=ii+1;
    
%     %% Check result of Three Channel FRET Computation against theory
%     % Only applicable for {Single Force, Single Unbinding Rate Constant, No Imaging Noise} Case
%     if strcmp(Name_DistributionType_F,'SingleValue') && strcmp(Name_DistributionType_kunbind,'SingleValue') && Toggle_Noise==0
%         F = F_base;
%         E0_Theory = MTS_FRETForceFunction(F); 
%         Eapp_Theory = E0_Theory./(1+num_D1A0./num_D1A1);
%         Sapp_Theory = (num_D1A1+num_D1A0)/(2*num_D1A1+num_D1A0+num_D0A1);
%         test_CompareSapp = abs(Sapp-Sapp_Theory);
%         test_CompareEapp = abs(Eapp-Eapp_Theory);
%         if (test_CompareSapp>1E-3) || (test_CompareEapp>1E-3)
%             warning( ['Result from Three Channel FRET Simulation disagrees with theory. [Sapp,Sapp_Theory,Eapp,Eapp_Theory] = [' num2str([Sapp,Sapp_Theory,Eapp,Eapp_Theory]) ']. UniqueID = ' UniqueID]);
%         end
%     end
  
    end
    disp(['Status: ' num2str(round(100*SimID/nSims,2)) '% complete.']);
    
    % Export SSAData 
    save(fullfile(AllPixelDataOutDir,['SSAData_' UniqueID '.mat']), 'SSAData');   

    % Export CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState as a table
    CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState = array2table(CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState,...
        'VariableNames',{'id_Pixel','Force','TurnoverRateConstant','State_Donor','State_Acceptor'});
    save(fullfile(AllPixelDataOutDir,['CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState_' UniqueID '.mat']), 'CombinedSensorPopulation_idPx_Force_UnbindRate_DState_AState');

    % Compute Stats for All Pixels for this SimID
    AllPixels_thisSimID = Data_OneRowPerPixel(Data_OneRowPerPixel.SimID==SimID,:);

    Force_mean_mean = nanmean(AllPixels_thisSimID.Force_mean);
    Lifetime_mean_mean = nanmean(AllPixels_thisSimID.Lifetime_mean);
    E0_mean_mean = nanmean(AllPixels_thisSimID.E0_mean);
    S0_mean_mean = nanmean(AllPixels_thisSimID.S0_mean);
    Eapp_mean = nanmean(AllPixels_thisSimID.Eapp);
    Sapp_mean = nanmean(AllPixels_thisSimID.Sapp);
    nD1A1_mean = nanmean(AllPixels_thisSimID.num_D1A1);
    nD1A0_mean = nanmean(AllPixels_thisSimID.num_D1A0);
    nD0A1_mean = nanmean(AllPixels_thisSimID.num_D0A1);
    nD0A0_mean = nanmean(AllPixels_thisSimID.num_D0A0);
    fractionDonorBroken_mean = nanmean(AllPixels_thisSimID.fractionDonorBroken);
    fractionAcceptorBroken_mean = nanmean(AllPixels_thisSimID.fractionAcceptorBroken);
    fractionSensorsBroken_mean = nanmean(AllPixels_thisSimID.fractionSensorsBroken);

    Force_mean_std = nanstd(AllPixels_thisSimID.Force_mean);
    Lifetime_mean_std = nanstd(AllPixels_thisSimID.Lifetime_mean);
    E0_mean_std = nanstd(AllPixels_thisSimID.E0_mean);
    S0_mean_std = nanstd(AllPixels_thisSimID.S0_mean);
    Eapp_std = nanstd(AllPixels_thisSimID.Eapp);
    Sapp_std = nanstd(AllPixels_thisSimID.Sapp);
    nD1A1_std = nanstd(AllPixels_thisSimID.num_D1A1);
    nD1A0_std = nanstd(AllPixels_thisSimID.num_D1A0);
    nD0A1_std = nanstd(AllPixels_thisSimID.num_D0A1);
    nD0A0_std = nanstd(AllPixels_thisSimID.num_D0A0);
    fractionDonorBroken_std = nanstd(AllPixels_thisSimID.fractionDonorBroken);
    fractionAcceptorBroken_std = nanstd(AllPixels_thisSimID.fractionAcceptorBroken);
    fractionSensorsBroken_std = nanstd(AllPixels_thisSimID.fractionSensorsBroken);

    Data_OneRowPerSimID(SimID,:) = table(...
        SimulationSetID,SimID,...
        Toggle_Noise,...
        LoadingRate,...
        BondType,...
        SensorUnbindingRateConstant_k01,SensorUnbindingRateConstant_xb1,...
        SensorUnbindingRateConstant_k02,SensorUnbindingRateConstant_xb2,...
        F_half_D,m_D,kms0_D,dxms_D,...
        F_half_A,m_A,kms0_A,dxms_A,...
        nSensors,nPixelsToSimulate,...
        Force_mean_mean,Lifetime_mean_mean,...
        E0_mean_mean,S0_mean_mean,Eapp_mean,Sapp_mean,nD1A1_mean,nD0A1_mean,nD1A0_mean,nD0A0_mean,fractionDonorBroken_mean,fractionAcceptorBroken_mean,fractionSensorsBroken_mean,...
        Force_mean_std,Lifetime_mean_std,...
        E0_mean_std,S0_mean_std,Eapp_std,Sapp_std,nD1A1_std,nD0A1_std,nD1A0_std,nD0A0_std,fractionDonorBroken_std,fractionAcceptorBroken_std,fractionSensorsBroken_std);
            
%     %% Check Mean of SSA against steady state solution
%     % Only applicable for {Single Force, Single Unbinding Rate Constant} Case
%     if strcmp(Name_DistributionType_F,'SingleValue') && strcmp(Name_DistributionType_kunbind,'SingleValue')
%         F = F_base;
%         kunload = SensorUnbindingRateConstant_k01.*exp(SensorUnbindingRateConstant_xb1.*F/kT) + ...
%             SensorUnbindingRateConstant_k02.*exp(SensorUnbindingRateConstant_xb2.*F/kT);
%         kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D,kT);
%         kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A,kT);
%         % Compute Steady State Prob of Each Species using Analytical Expressions
%         nD1A1_Steady = nSensors .* kunload/(kmsA + kmsD + kunload);
%         nD0A1_Steady = nSensors .* (kmsD*kunload)/((kmsA + kunload)*(kmsA + kmsD + kunload));
%         nD1A0_Steady = nSensors .* (kmsA*kunload)/((kmsD + kunload)*(kmsA + kmsD + kunload));
%         nD0A0_Steady = nSensors .* (kmsA*kmsD*(kmsA + kmsD + 2*kunload))/((kmsA + kunload)*(kmsD + kunload)*(kmsA + kmsD + kunload));
%         
%         test_CompareToSteady_nD1A1 = abs(nD1A1_mean-nD1A1_Steady);
%         test_CompareToSteady_nD0A1 = abs(nD0A1_mean-nD0A1_Steady);
%         test_CompareToSteady_nD1A0 = abs(nD1A0_mean-nD1A0_Steady);
%         test_CompareToSteady_nD0A0 = abs(nD0A0_mean-nD0A0_Steady);
%         test_CompareToSteady_max = max([test_CompareToSteady_nD1A1,test_CompareToSteady_nD0A1,test_CompareToSteady_nD1A0,test_CompareToSteady_nD0A0]);
%         % test_CompareToSteady = sum(abs([nD1A1_mean,nD0A1_mean,nD1A0_mean,nD0A0_mean]-[nD1A1_Steady,nD0A1_Steady,nD1A0_Steady,nD0A0_Steady]));
%         if test_CompareToSteady_max > 1
%             warning( ['Result from SSA for [nD1A1_mean,nD0A1_mean,nD1A0_mean,nD0A0_mean] does not match steady state solution. UniqueID = ' UniqueID '. test_CompareToSteady = ' num2str([test_CompareToSteady_nD1A1,test_CompareToSteady_nD0A1,test_CompareToSteady_nD1A0,test_CompareToSteady_nD0A0])] );
%         end
%     end
        
end
end
end

save(fullfile(AllPixelDataOutDir,['Data_OneRowPerSimID_' StudyName]), 'Data_OneRowPerSimID');
% save(fullfile(PlotOutDir,['Data_OneRowPerSimID_' StudyName]), 'Data_OneRowPerSimID');

save(fullfile(AllPixelDataOutDir,['Data_OneRowPerPixel_' StudyName]), 'Data_OneRowPerPixel');
% save(fullfile(PlotOutDir,['Data_OneRowPerPixel_' StudyName]), 'Data_OneRowPerPixel');

%% Figures
load('../../TheoreticalESHistogramsForMTS\ESCurves_TSMod_Tension.mat');
addpath('PlotFunctions_ESHistograms');
set(0,'defaultAxesFontSize',10);

%% Settings for {Eapp,Sapp} Histrograms
% Bin Edges for Fine ES-Histogram
% FixedBinEdges_Eff = [0:.005:.5];
% FixedBinEdges_Stoich = [0:.010:1];
FixedBinEdges_Eff = [0:.01:.5];
FixedBinEdges_Stoich = [0:.020:1];

% Bin Edges for Coarse ES-Histogram
FixedBinEdges_Eff_Coarse = [0:.025:0.5];
FixedBinEdges_Stoich_Coarse = [0:.05:1];

% Handling of Out-of-Bounds E and S
Method_HandleOutOfBoundsEandS = 'SetToLimit'; % 'SetToLimit' = set to min or max, 'Exclude' = exclude both E and S

% Min/Max
eff_min = 0;
eff_max = .5;
stoic_min = 0;
stoic_max = 1;

%% Settings for {Idd/Iaa,Fc/Iaa} Histograms
% Bin Edges
FixedBinEdges_Idd_per_Iaa = linspace(0,3.5,100);
FixedBinEdges_Fc_per_Iaa = linspace(0,3.5,100);

% Handling of Out-of-Bounds E and S
Method_HandleOutOfBoundsIddperIaaAndFcperIaa = 'SetToLimit'; % 'SetToLimit' = set to min or max, 'Exclude' = exclude both E and S

% Min/Max
Idd_per_Iaa_min = 0;
Idd_per_Iaa_max = 3.5;
Fc_per_Iaa_min = 0;
Fc_per_Iaa_max = 3.5;

%% Make ES Histograms
ShowEmptyBins = 'off';
Toggle_GridLines = 0;
for SimID = [unique(Data_OneRowPerSimID.SimID)']
   
    % Load Px Data
    AllPixels_thisSimID = Data_OneRowPerPixel(Data_OneRowPerPixel.SimID==SimID,:);
    Eapp = AllPixels_thisSimID.Eapp;
    Sapp = AllPixels_thisSimID.Sapp;
    Idd_per_Iaa = AllPixels_thisSimID.Idd./AllPixels_thisSimID.Iaa;
    Fc_per_Iaa = AllPixels_thisSimID.Fc./AllPixels_thisSimID.Iaa;
    
    % Apply Method to Handle Out of Bounds Eapp,Sapp
    if strcmp(Method_HandleOutOfBoundsEandS,'SetToLimit')
        Eapp(Eapp<eff_min) = eff_min;
        Eapp(Eapp>eff_max) = eff_max;
        Sapp(Sapp<stoic_min) = stoic_min;
        Sapp(Sapp>stoic_max) = stoic_max;
    elseif strcmp(Method_HandleOutOfBoundsEandS,'Exclude')
        Exclude = Eapp<eff_min | Eapp>eff_max | Sapp<stoic_min | Sapp>stoic_max;
        Eapp(Exclude) = [];
        Sapp(Exclude) = [];
    end
    
    % Apply Method to Handle Out of Bounds Idd/Iaa,Fc/Iaa
    if strcmp(Method_HandleOutOfBoundsIddperIaaAndFcperIaa,'SetToLimit')
        Idd_per_Iaa(Idd_per_Iaa<Idd_per_Iaa_min) = Idd_per_Iaa_min;
        Idd_per_Iaa(Idd_per_Iaa>Idd_per_Iaa_max) = Idd_per_Iaa_max;
        Fc_per_Iaa(Fc_per_Iaa<Fc_per_Iaa_min) = Fc_per_Iaa_min;
        Fc_per_Iaa(Fc_per_Iaa>Fc_per_Iaa_max) = Fc_per_Iaa_max;
    elseif strcmp(Method_HandleOutOfBoundsIddperIaaAndFcperIaa,'Exclude')
        Exclude = Idd_per_Iaa<Idd_per_Iaa_min | Idd_per_Iaa>Idd_per_Iaa_max | Fc_per_Iaa<Fc_per_Iaa_min | Fc_per_Iaa>Fc_per_Iaa_max;
        Idd_per_Iaa(Exclude) = [];
        Fc_per_Iaa(Exclude) = [];
    end

    % Make ES Histogram
    f=figure(SimID);
    % Set Eapp bins (x variable)
    edgesX = FixedBinEdges_Eff;
    % Set Sapp bins (y variable)
    edgesY = FixedBinEdges_Stoich;
    % Plot 2D histogram using MATLAB function
    h = histogram2(Eapp,Sapp,edgesX,edgesY,'DisplayStyle','tile','ShowEmptyBins',ShowEmptyBins); hold on;
    colormap('parula');
    xlabel('E');
    ylabel('S');
    xlim([0 .5]);
    ylim([0 1]);
    set(gca,'xtick',[0:.05:.5]);
    set(gca,'ytick',[0:.1:1]);
    c = colorbar('Location','south');
    if Toggle_GridLines
        ax = gca;
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = FixedBinEdges_Eff;
        ax.YAxis.MinorTick = 'on';
        ax.YAxis.MinorTickValues = FixedBinEdges_Stoich;
        set(gca, 'xminorgrid', 'on')
        set(gca, 'yminorgrid', 'on')
        set(gca,'MinorGridLineStyle','-')
    end
    
    % Plot ES Curves
    if strcmp(ShowEmptyBins,'off')
        RefLineStyle = 'k--';
    else
        RefLineStyle = 'w--';
    end
    for ff=[1 4 7] % ff=[1 4 7 10 13 16] for F=0,3,6,9,12,15
        p = plot(ESCurves{ff}.Eapp,ESCurves{ff}.Sapp,RefLineStyle);
        hold on;
    end
    
    % S0 Line
    plot([0 .5], [ESCurves{1}.S0(1) ESCurves{1}.S0(1)], RefLineStyle);
    
    % Unloaded E0 Line
    plot([ESCurves{1}.E0(1) ESCurves{1}.E0(1)], [0 1], RefLineStyle);
    
    % Save full figure panel in specified file formats
    FigName = ['EShist_SimID_' num2str(SimID)];
    Resolution_png = 300;
    print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
    print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');

    close(f);
end

%% Plot Means +/- Std Dev and Compare to Steady Solution **NEW**   
load('C:\Users\tcsho\OneDrive - Duke University\01 Code Repository\03 Models\FPMechSwitch_v02/TheoreticalESHistogramsForMTS\ESCurves_TSMod_Tension.mat');
set(0,'defaultAxesFontSize',16);

SimIDSets = {{0*7+[1:7]},{1*7+[1:7]},{2*7+[1:7]}; % Ideal/Slip/Catch-slip for No FP Break
    {3*7+[1:7]},{4*7+[1:7]},{5*7+[1:7]}; % Ideal/Slip/Catch-slip for A Break
    {6*7+[1:7]},{7*7+[1:7]},{8*7+[1:7]}; % Ideal/Slip/Catch-slip for D Break
    {9*7+[1:7]},{10*7+[1:7]},{11*7+[1:7]}}; % Ideal/Slip/Catch-slip Both Break

FigSuffix = {'None','AccOnly', 'DonOnly', 'Both'};

figID = 0;
for row = 1:length(SimIDSets(:,1))
    figID=figID+1;
    f = figure;
for ss=1:length(SimIDSets(row,:))
    
    SimIDs = SimIDSets{row,ss}; SimIDs = SimIDs{:};

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

    % Get Mean and Std Dev for Eapp and Sapp from Stochastic Simulations
    thisData = Data_OneRowPerSimID(SimIDs,:);
    Eapp_mean = thisData.Eapp_mean;
    Sapp_mean = thisData.Sapp_mean;
    Eapp_std = thisData.Eapp_std;
    Sapp_std = thisData.Sapp_std;
    
    % Plot ErrorBar of Sim Data
    LineStyle = '.-';
    MarkerSize = 15;
    LineWidth = 1.5;
    % errorbar(x,y,yneg,ypos,xneg,xpos)
    eb = errorbar(Eapp_mean,Sapp_mean,Sapp_std,Sapp_std,Eapp_std,Eapp_std,LineStyle,'MarkerSize',MarkerSize);
    eb.LineWidth = LineWidth;
    hold on;
    
end
    FigName = ['ES_SSAMeansStdDevs_LoadRate_ThreeBondTypes_' FigSuffix{figID}];
    Resolution_png = 300;
    print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
    print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');
    % close(f);
end

%% Plot Means and Compare to Steady Solution **NEW**   
load('C:\Users\tcsho\OneDrive - Duke University\01 Code Repository\03 Models\FPMechSwitch_v02/TheoreticalESHistogramsForMTS\ESCurves_TSMod_Tension.mat');
set(0,'defaultAxesFontSize',16);

SimIDSets = {{0*7+[1:7]},{1*7+[1:7]},{2*7+[1:7]}; % Ideal/Slip/Catch-slip for No FP Break
    {3*7+[1:7]},{4*7+[1:7]},{5*7+[1:7]}; % Ideal/Slip/Catch-slip for A Break
    {6*7+[1:7]},{7*7+[1:7]},{8*7+[1:7]}; % Ideal/Slip/Catch-slip for D Break
    {9*7+[1:7]},{10*7+[1:7]},{11*7+[1:7]}}; % Ideal/Slip/Catch-slip Both Break

FigSuffix = {'None','AccOnly', 'DonOnly', 'Both'};

figID = 0;
for row = 1:length(SimIDSets(:,1))
    figID=figID+1;
    f = figure;
for ss=1:length(SimIDSets(row,:))
    
    SimIDs = SimIDSets{row,ss}; SimIDs = SimIDs{:};

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

    % Get Mean and Std Dev for Eapp and Sapp from Stochastic Simulations
    thisData = Data_OneRowPerSimID(SimIDs,:);
    Eapp_mean = thisData.Eapp_mean;
    Sapp_mean = thisData.Sapp_mean;
    Eapp_std = thisData.Eapp_std;
    Sapp_std = thisData.Sapp_std;
    
    % Plot ErrorBar of Sim Data
    LineStyle = '.-';
    MarkerSize = 15;
    LineWidth = 1.5;
    plot(Eapp_mean,Sapp_mean,LineStyle,'MarkerSize',MarkerSize);
    hold on;
    
end
    FigName = ['ES_SSAMeans_LoadRate_ThreeBondTypes_' FigSuffix{figID}];
    Resolution_png = 300;
    print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
    print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');
    % close(f);
end

%% Support Functions
% FP Mechanial Switching, Forward Rate Constant
function out = kms(F,F_half,m,kms0,dxms,kT)
    out = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
end