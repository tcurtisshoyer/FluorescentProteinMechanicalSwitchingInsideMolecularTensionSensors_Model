close all; clear all; clc;

% Goal: Validate derived expressions for Eapp and Sapp in terms of distribution
% of sensor states [nD1A1,nD0A1,nD1A0,nD0A0] by comparison to Eapp and Sapp
% computed directly from channel intensities (see
% simulateThreeChannelFRET_E0Dist.m). The derived expressions hold for the
% case of a single force F, and thus a single E0=f(F). Also, nD1A1 must be
% > 0.

%% Derived Expressions for Eapp and Sapp to Validate
% Eapp = E0./(1+nD1A0./nD1A1);
% Sapp = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);

%% Explicit Computation of Eapp and Sapp from Channel Intensities
% See: simulateThreeChannelFRET_E0Dist.m
% Ground truth calibration constants
% G_true = 1.65;
% k_true = 1.00;
% gammaM_true = G_true;
% betaX_true = 1/G_true/k_true;
% 
% Donor Bleedthru and Acceptor Direct Excitation constants
% alphaBT_true = @(x) 0.75; % Actual donor bleedthru coefficient. Can be constant or function of intensity, x=Idd_TOT.
% deltaDE_true = @(x) 0.25; % Actual acceptor direct excitation coefficient. Can be constant or function of intensity, x=Idd_TOT.
% alphaBT_est = 0.75; % Estimated donor bleedthru coefficient using donor only experiments and intensity-independent assumption
% deltaDE_est = 0.25; % Estimated acceptor direct excitation coefficient using acceptor only experiments and intensity-independent assumption
% 
% Compute Total Channel Intensities for the general mixed population of
% sensor states based on numbers of sensors in each state, ideal FRET Eff
% of each sensor, and and the ground truth calibration constants.
% k = 100; % mean signal in AA channel per acceptor fluorophore
% Iaa_TOT = sum( (state_D1A1+state_D0A1) .* k );
% Idd_TOT = sum( (state_D1A1.*(1-E0) + state_D1A0) .* k/gammaM_true/betaX_true );
% Ida_TOT = sum( E0.*state_D1A1 .*k/betaX_true ) + alphaBT_true(Idd_TOT).*Idd_TOT + deltaDE_true(Iaa_TOT).*Iaa_TOT;
%         
% Compute corrected FRET intensity (Fc or Ida_corr) using estimated
% coefficients for donor bleedthru and acceptor direct excitation
% Fc = Ida - alphaBT_est*Idd - deltaDE_est*Iaa;
% 
% Compute Apparent FRET Eff (Eapp)
% Eapp = Fc./(gammaM_true*Idd + Fc);
% 
% Comptue Apparent Stoichiometry (Sapp)
% Sapp = (gammaM_true*Idd + Fc)./(gammaM_true*Idd + Fc + Iaa/betaX_true);

%% MTS FRET Eff vs Force Model
% Time-Averaged Signal of Single MTS at Specified Force
load('../TheoreticalESHistogramsForMTS\LUT_OriginalTSMod.mat');
MTS_FRETForceFunction = @(F) max(interp1(LUT_OriginalTSMod.Force_app, LUT_OriginalTSMod.FRETEff, F),0);

%% Run Validation - Sweep nD1A0
StudyID = 1;
AllData{StudyID} = [];
nSensors = 50;
for F=[0 3 6 9] 
for nD1A0 = [0:3:48]
for nD0A1 = [0]
for nD0A0 = [0]
    
    nD1A1 = nSensors - (nD0A1+nD1A0+nD0A0);
    
    if nD1A1<0
        error('Sensor State Distribution invalid.');
    end
    
    %% Compute Eapp and Sapp from Derived Expressions
    % Use analytical expression that requires all sensors in ensemble have the
    % same force and assumes that D1A1>0
    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 
    % Compute Eapp
    Eapp_Derived = E0./(1+nD1A0./nD1A1);
    % Compute Sapp
    S0 = 0.5;
    Sapp_Derived = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);
    
    %% Compute Eapp and Sapp Directly From Channel Intensities

    % Initiate the Sensor Population with Single Force
    SensorForceDistribution = F.*ones(nSensors,1);
    AllSensors_State_Donor = [];
    AllSensors_State_Acceptor = [];
    for ii=1:nD1A1
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+1:nD1A1+nD1A0
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 0;
    end
    for ii=nD1A1+nD1A0+1:nD1A1+nD1A0+nD0A1
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+nD1A0+nD0A1+1:nSensors
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 0;
    end

    % Compute Three Channel FRET Measurements of MTS Ensemble            
    % Obtain FRET Eff for Ideal Sensor State for Each Sensor based
    % on Sensor Force Distribution
    AllSensors_E0 = MTS_FRETForceFunction(SensorForceDistribution);
    AllSensors_S0 = 0.5*ones(size(AllSensors_E0));
    % Store D#A# State for Each Sensor to pass to
    % simulateThreeChannelFRET_E0dist.m
    AllSensors_D1A1 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1; % Set to 1 if D1A1, else 0
    AllSensors_D1A0 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0; % Set to 1 if D1A0, else 0
    AllSensors_D0A1 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1; % Set to 1 if D0A1, else 0
    AllSensors_D0A0 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0; % Set to 1 if D0A0, else 0
    % Run Computations
    Toggle_Noise = 0;
    PixelData = simulateThreeChannelFRET_E0dist(AllSensors_D1A1, AllSensors_D1A0, AllSensors_D0A1, AllSensors_D0A0, AllSensors_E0, Toggle_Noise);
    Iaa = PixelData.Iaa; 
    Idd = PixelData.Idd;
    Ida = PixelData.Ida;
    Fc = PixelData.Fc;
    Eapp = PixelData.Eapp;
    Sapp = PixelData.Sapp;
    
    %% Log Data
    AllData{StudyID} = [AllData{StudyID}; table(nD1A1,nD1A0,nD0A1,nD0A0,F,E0,Eapp_Derived,Sapp_Derived,Iaa,Idd,Ida,Fc,Eapp,Sapp)];
    
end
end
end
end

%% Run Validation - Sweep nD0A1
StudyID = 2;
AllData{StudyID} = [];
nSensors = 50;
for F=[0 3 6 9] 
for nD1A0 = [0]
for nD0A1 = [0:3:48]
for nD0A0 = [0]
    
    nD1A1 = nSensors - (nD0A1+nD1A0+nD0A0);
    
    if nD1A1<0
        error('Sensor State Distribution invalid.');
    end
    
    %% Compute Eapp and Sapp from Derived Expressions
    % Use analytical expression that requires all sensors in ensemble have the
    % same force and assumes that D1A1>0
    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 
    % Compute Eapp
    Eapp_Derived = E0./(1+nD1A0./nD1A1);
    % Compute Sapp
    S0 = 0.5;
    Sapp_Derived = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);
    
    %% Compute Eapp and Sapp Directly From Channel Intensities

    % Initiate the Sensor Population with Single Force
    SensorForceDistribution = F.*ones(nSensors,1);
    AllSensors_State_Donor = [];
    AllSensors_State_Acceptor = [];
    for ii=1:nD1A1
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+1:nD1A1+nD1A0
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 0;
    end
    for ii=nD1A1+nD1A0+1:nD1A1+nD1A0+nD0A1
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+nD1A0+nD0A1+1:nSensors
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 0;
    end

    % Compute Three Channel FRET Measurements of MTS Ensemble            
    % Obtain FRET Eff for Ideal Sensor State for Each Sensor based
    % on Sensor Force Distribution
    AllSensors_E0 = MTS_FRETForceFunction(SensorForceDistribution);
    AllSensors_S0 = 0.5*ones(size(AllSensors_E0));
    % Store D#A# State for Each Sensor to pass to
    % simulateThreeChannelFRET_E0dist.m
    AllSensors_D1A1 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1; % Set to 1 if D1A1, else 0
    AllSensors_D1A0 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0; % Set to 1 if D1A0, else 0
    AllSensors_D0A1 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1; % Set to 1 if D0A1, else 0
    AllSensors_D0A0 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0; % Set to 1 if D0A0, else 0
    % Run Computations
    Toggle_Noise = 0;
    PixelData = simulateThreeChannelFRET_E0dist(AllSensors_D1A1, AllSensors_D1A0, AllSensors_D0A1, AllSensors_D0A0, AllSensors_E0, Toggle_Noise);
    Iaa = PixelData.Iaa; 
    Idd = PixelData.Idd;
    Ida = PixelData.Ida;
    Fc = PixelData.Fc;
    Eapp = PixelData.Eapp;
    Sapp = PixelData.Sapp;
    
    %% Log Data
    AllData{StudyID} = [AllData{StudyID}; table(nD1A1,nD1A0,nD0A1,nD0A0,F,E0,Eapp_Derived,Sapp_Derived,Iaa,Idd,Ida,Fc,Eapp,Sapp)];
    
end
end
end
end

%% Run Validation - Sweep nD0A0
StudyID = 3;
AllData{StudyID} = [];
nSensors = 50;
for F=[0 3 6 9] 
for nD1A0 = [0]
for nD0A1 = [0]
for nD0A0 = [0:3:48]
    
    nD1A1 = nSensors - (nD0A1+nD1A0+nD0A0);
    
    if nD1A1<0
        error('Sensor State Distribution invalid.');
    end
    
    %% Compute Eapp and Sapp from Derived Expressions
    % Use analytical expression that requires all sensors in ensemble have the
    % same force and assumes that D1A1>0
    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 
    % Compute Eapp
    Eapp_Derived = E0./(1+nD1A0./nD1A1);
    % Compute Sapp
    S0 = 0.5;
    Sapp_Derived = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);
    
    %% Compute Eapp and Sapp Directly From Channel Intensities

    % Initiate the Sensor Population with Single Force
    SensorForceDistribution = F.*ones(nSensors,1);
    AllSensors_State_Donor = [];
    AllSensors_State_Acceptor = [];
    for ii=1:nD1A1
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+1:nD1A1+nD1A0
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 0;
    end
    for ii=nD1A1+nD1A0+1:nD1A1+nD1A0+nD0A1
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+nD1A0+nD0A1+1:nSensors
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 0;
    end

    % Compute Three Channel FRET Measurements of MTS Ensemble            
    % Obtain FRET Eff for Ideal Sensor State for Each Sensor based
    % on Sensor Force Distribution
    AllSensors_E0 = MTS_FRETForceFunction(SensorForceDistribution);
    AllSensors_S0 = 0.5*ones(size(AllSensors_E0));
    % Store D#A# State for Each Sensor to pass to
    % simulateThreeChannelFRET_E0dist.m
    AllSensors_D1A1 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1; % Set to 1 if D1A1, else 0
    AllSensors_D1A0 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0; % Set to 1 if D1A0, else 0
    AllSensors_D0A1 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1; % Set to 1 if D0A1, else 0
    AllSensors_D0A0 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0; % Set to 1 if D0A0, else 0
    % Run Computations
    Toggle_Noise = 0;
    PixelData = simulateThreeChannelFRET_E0dist(AllSensors_D1A1, AllSensors_D1A0, AllSensors_D0A1, AllSensors_D0A0, AllSensors_E0, Toggle_Noise);
    Iaa = PixelData.Iaa; 
    Idd = PixelData.Idd;
    Ida = PixelData.Ida;
    Fc = PixelData.Fc;
    Eapp = PixelData.Eapp;
    Sapp = PixelData.Sapp;
    
    %% Log Data
    AllData{StudyID} = [AllData{StudyID}; table(nD1A1,nD1A0,nD0A1,nD0A0,F,E0,Eapp_Derived,Sapp_Derived,Iaa,Idd,Ida,Fc,Eapp,Sapp)];
    
end
end
end
end

%% Run Validation - Sweep nD1A0=nD0A1
StudyID = 4;
AllData{StudyID} = [];
nSensors = 50;
for F=[0 3 6 9] 
for nD1A0 = [0:2:24]
for nD0A1 = nD1A0
for nD0A0 = [0]
    
    nD1A1 = nSensors - (nD0A1+nD1A0+nD0A0);
    
    if nD1A1<0
        error('Sensor State Distribution invalid.');
    end
    
    %% Compute Eapp and Sapp from Derived Expressions
    % Use analytical expression that requires all sensors in ensemble have the
    % same force and assumes that D1A1>0
    % Obtain FRET Eff for Ideal Sensor State based on FRET-Force Calibration for MTS
    E0 = MTS_FRETForceFunction(F); 
    % Compute Eapp
    Eapp_Derived = E0./(1+nD1A0./nD1A1);
    % Compute Sapp
    S0 = 0.5;
    Sapp_Derived = (nD1A1+nD1A0)/(2*nD1A1+nD1A0+nD0A1);
    
    %% Compute Eapp and Sapp Directly From Channel Intensities

    % Initiate the Sensor Population with Single Force
    SensorForceDistribution = F.*ones(nSensors,1);
    AllSensors_State_Donor = [];
    AllSensors_State_Acceptor = [];
    for ii=1:nD1A1
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+1:nD1A1+nD1A0
        AllSensors_State_Donor(ii,:) = 1;
        AllSensors_State_Acceptor(ii,:) = 0;
    end
    for ii=nD1A1+nD1A0+1:nD1A1+nD1A0+nD0A1
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 1;
    end
    for ii=nD1A1+nD1A0+nD0A1+1:nSensors
        AllSensors_State_Donor(ii,:) = 0;
        AllSensors_State_Acceptor(ii,:) = 0;
    end

    % Compute Three Channel FRET Measurements of MTS Ensemble            
    % Obtain FRET Eff for Ideal Sensor State for Each Sensor based
    % on Sensor Force Distribution
    AllSensors_E0 = MTS_FRETForceFunction(SensorForceDistribution);
    AllSensors_S0 = 0.5*ones(size(AllSensors_E0));
    % Store D#A# State for Each Sensor to pass to
    % simulateThreeChannelFRET_E0dist.m
    AllSensors_D1A1 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1; % Set to 1 if D1A1, else 0
    AllSensors_D1A0 = AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0; % Set to 1 if D1A0, else 0
    AllSensors_D0A1 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1; % Set to 1 if D0A1, else 0
    AllSensors_D0A0 = AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0; % Set to 1 if D0A0, else 0
    % Run Computations
    Toggle_Noise = 0;
    PixelData = simulateThreeChannelFRET_E0dist(AllSensors_D1A1, AllSensors_D1A0, AllSensors_D0A1, AllSensors_D0A0, AllSensors_E0, Toggle_Noise);
    Iaa = PixelData.Iaa; 
    Idd = PixelData.Idd;
    Ida = PixelData.Ida;
    Fc = PixelData.Fc;
    Eapp = PixelData.Eapp;
    Sapp = PixelData.Sapp;
    
    %% Log Data
    AllData{StudyID} = [AllData{StudyID}; table(nD1A1,nD1A0,nD0A1,nD0A0,F,E0,Eapp_Derived,Sapp_Derived,Iaa,Idd,Ida,Fc,Eapp,Sapp)];
    
end
end
end
end

%% Save All Data
PlotOutDir = 'Validate_EappSappDerivedExpressions';
save(fullfile(PlotOutDir,'AllData.mat'), 'AllData');

%% Plot Data
xVars = {'nD1A0','nD0A1','nD0A0','nD1A0'};
for StudyID = 1:4
    xVar = xVars{StudyID};
    f = figure;
    f.Position = [50 50 300 700];
    
    subplot(3,1,2);
    MarkerSize = [4 4 4 4];
    ii=1;
    for F = [0 3 6 9]
        thisData = AllData{StudyID}(AllData{StudyID}.F==F,:);
        p = plot(thisData.(xVar),thisData.Eapp,'o','MarkerSize',MarkerSize(ii));
        hold on;
        plot(thisData.(xVar),thisData.Eapp_Derived,'-','LineWidth',1,'Color',p.Color);
        ii=ii+1;
    end
    ylabel('Eapp');
    xlabel(xVar);
    set(gca, 'FontName', 'Arial');
    subplot(3,1,3);
    MarkerSize = [6 5 4 3];
    ii=1;
    for F = [0 3 6 9]
        thisData = AllData{StudyID}(AllData{StudyID}.F==F,:);
        p = plot(thisData.(xVar),thisData.Sapp,'o','MarkerSize',MarkerSize(ii));
        hold on;
        plot(thisData.(xVar),thisData.Sapp_Derived,'-','LineWidth',1,'Color','k');
        ii=ii+1;
    end
    ylabel('Sapp'); ylim([0 1]);
    xlabel(xVar);
    set(gca, 'FontName', 'Arial');
    subplot(3,1,1);
    ii=1;
    for F = [3]
        thisData = AllData{StudyID}(AllData{StudyID}.F==F,:);
        plot(thisData.(xVar),thisData.Iaa,'o','MarkerSize',4,'Color',"#EDB120");
        hold on;
        plot(thisData.(xVar),thisData.Idd,'o','MarkerSize',4,'Color',"#0072BD");
        plot(thisData.(xVar),thisData.Fc,'o','MarkerSize',4,'Color',"#D95319");
        ii=ii+1;
    end
    ylabel('I'); ylim([0 6000]);
    xlabel(xVar);
    % legend('Iaa, F=3','Idd, F=3','Fc, F=3','Location','Best');
    set(gca, 'FontName', 'Arial');
    
    FigName = ['Plots_StudyID_' num2str(StudyID)];
    Resolution_png = 300;
    print(f,fullfile(PlotOutDir,[FigName '.png']),'-dpng',['-r' num2str(Resolution_png)],'-painters');
    print(f,fullfile(PlotOutDir,[FigName '.svg']),'-dsvg','-painters');
    
end