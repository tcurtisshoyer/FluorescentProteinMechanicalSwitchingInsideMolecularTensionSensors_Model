function [AllSensors_ActualForce,AllSensors_Lifetime,AllSensors_State_Donor,AllSensors_State_Acceptor,num_D1A1,num_D1A0,num_D0A1,num_D0A0,SimulationPerformanceStats,DataLog_Ensemble,DataLog_Transitions,DataLog_LifetimesAtTurnoverOnly,AllVariableNames] = SSA_4StateSensor_ForceDistribution_ForceDepBonds(ModelParams,SimulationParams,Toggle_OutputTimeCourseOfEnsemble,Toggle_OutputLogOfAllTransitions)

%% Simulation Parameters
tend = SimulationParams.tend; % Total simulation time
interval_LogEnsembleData = SimulationParams.interval_LogEnsembleData; % Time interval for logging linkage data 
timeSpan_LogTransitionData = SimulationParams.timeSpan_LogTransitionData; % Time span over which to log linkage transisitions
maxAllowedTimeStep = SimulationParams.maxAllowedTimeStep; % Max Allowed Timestep in Gillespie SSA. If dt larger than this value, advance time and update mechanics without executing a linkage transition
rateMax = 1E15; % Max value for force-dependent rate constant evaluation
kT=4.114; % [pN/nm] Approximate value of kT at 298 K

%% Model Parameters

% Number of Sensors
numSensors = ModelParams.numSensors;

% Sensor Force Distribution - Each "sensor binding site" has specified force
SensorForceDistribution = ModelParams.SensorForceDistribution;

% Mech Switching Paramters for Donor - All sensors have same kinetics
% kmsD = kms(F,F_half_D,m_D,kms0_D,dxms_D);
F_half_D = ModelParams.F_half_D;
m_D = ModelParams.m_D;
kms0_D = ModelParams.kms0_D;
dxms_D = ModelParams.dxms_D;

% Mech Switching Paramters for Acceptor - All sensors have same kinetics
% kmsA = kms(F,F_half_A,m_A,kms0_A,dxms_A);
F_half_A = ModelParams.F_half_A;
m_A = ModelParams.m_A;
kms0_A = ModelParams.kms0_A;
dxms_A = ModelParams.dxms_A;

% Sensor Unbinding Rate Constant - Specified as Parameters for 2 Pathway
% Catch-Slip Bond Model: koff(F)k01.*exp(x1.*F/kT) + k02.*exp(x2.*F/kT)
SensorUnbindingRateConstant_k01 = ModelParams.SensorUnbindingRateConstant_k01;
SensorUnbindingRateConstant_xb1 = ModelParams.SensorUnbindingRateConstant_xb1;
SensorUnbindingRateConstant_k02 = ModelParams.SensorUnbindingRateConstant_k02;
SensorUnbindingRateConstant_xb2 = ModelParams.SensorUnbindingRateConstant_xb2;

%% Initial Conditions
AllSensors_ActualForce = SensorForceDistribution;
AllSensors_State_Donor = ones(numSensors,1); % FP State: 1=On, 0=Off
AllSensors_State_Acceptor = ones(numSensors,1); % FP State: 1=On, 0=Off

%% Structure to track time since last turnover for each sensor
AllSensors_Lifetime = zeros(numSensors,1);

%% Initialize Time
t = 0;
dt = 0;

%% Pre-Allocate Data Storage
% DataLog_Transitions: A data structure to log all transition events. Logs the time of occurence, id for
% the linkage that underwent the transition, id of the transition, and
% force across the sensor.
initHeight = 1E5;
VariableNames_DataLog_Transitions = {'t','idSensor','idStateVariable','StateValueFrom','StateValueTo','RateConstant','RateConstantForce','GillespieActivity','TimeSinceLastTurnover'};
PreAllocationBlocks.DataLog_Transitions = zeros(initHeight,length(VariableNames_DataLog_Transitions));
DataLog_Transitions = PreAllocationBlocks.DataLog_Transitions;
iDataLog_Transitions = 0; % Index of current data log. Incremented in simulation.

% DataLog_Ensemble: A data structure to log time course data for a
% the whole sensor ensemble
initHeight = 1E5;
VariableNames_DataLog_Ensemble = {'t',...
        'numD1A1',...
        'numD1A0',...
        'numD0A1',...
        'numD0A0'};
PreAllocationBlocks.DataLog_Ensemble = zeros(initHeight,length(VariableNames_DataLog_Ensemble));
DataLog_Ensemble = PreAllocationBlocks.DataLog_Ensemble;
iDataLog_Ensemble = 0; % Index of current data log. Incremented in simulation.
t_lastDataLogEnsemble = t;

%% Stochastic Simulation Algorithm
numAcceptedTransitions = 0; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt, so a reaction will occur along with update to forces
numRejectedTransitions = 0; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt, so no reaction will occur, just update forces
avgTimeStep = 0; % Record simulation performance stats: avg length of time step
DataLog_LifetimesAtTurnoverOnly = []; % Stores the lifetime of linkages when they have a turnover event
while t<tend
    
    %% [1] Execute Next Transition and Update Time (Gillespie Direct Method)
    % Summary of Approach:
    % [1A] Calculate rate constants for each available FP state transitions in each MTS based on the kinetic models for FPs and the force across the MTS.
    % [1B] Determine and execute earliest state transition and update time
    % using the Direct Gillespie SSA.

    % [1A] Find all available transitions (all transitions open to each
    % Sensor) and compute each transition rate.
    
    % Possible Transitions is a Nx6 matrix where each row is a possible
    % transition available at the current point in time:
    % Col 1: Sensor ID
    % Col 2: State Variable ID:
        % State_Donor = 1
        % State_Acceptor = 2
        % NOTE: Turnover Transition (D#A#-->D1A1) indicated using State ID of 0 and State Value 0-->0
    % Col 3: Current/From State Variable Value
        % NOTE: Turnover Transition (D#A#-->D1A1) indicated using State ID of 0 and State Value 0-->0
    % Col 4: New/To State Variable Value
        % NOTE: Turnover Transition (D#A#-->D1A1) indicated using State ID of 0 and State Value 0-->0
    % Col 5: Rate Constant
    % Col 6: Force used to compute Rate Constant
    % 
    % Example for acceptor mechanical switching forward transition:
    % 1|2|1|0|kms(F,F_half_A,m_A,kms0_A,dxms_A)|F

    allTransitions = zeros(numSensors*10,6); % Initialize allTransitions to save time
    count=0;
    for thisSensorID = 1:numSensors
        F = AllSensors_ActualForce(thisSensorID);
        State_Donor = AllSensors_State_Donor(thisSensorID);
        State_Acceptor = AllSensors_State_Acceptor(thisSensorID);

        % Possible Transitions for Donor FP
        thisStateVariableID = 1; % Donor FP
        currentStateValue = State_Donor;
        if currentStateValue==1 % On/Folded --> Off/Unfolded
            rate = kms(F,F_half_D,m_D,kms0_D,dxms_D);
            thisPossibleTransition = [thisSensorID thisStateVariableID 1 0 rate F];
            allTransitions(count+1,:) = thisPossibleTransition;
            count = count + 1;
        end
        
        % Possible Transitions for Acceptor FP
        thisStateVariableID = 2; % Acceptor FP
        currentStateValue = State_Acceptor;
        if currentStateValue==1 % On/Folded --> Off/Unfolded
            rate = kms(F,F_half_A,m_A,kms0_A,dxms_A);
            thisPossibleTransition = [thisSensorID thisStateVariableID 1 0 rate F];
            allTransitions(count+1,:) = thisPossibleTransition;
            count = count + 1;
        end
        
        % Turnover Transition
        % All sensors subject to turnover
        % NOTE: Turnover Transition (D#A#-->D1A1) indicated using State ID of 0 and State Value 0-->0
        thisStateVariableID = 0;
        koff = SensorUnbindingRateConstant_k01.*exp(SensorUnbindingRateConstant_xb1.*F/kT) + ...
            SensorUnbindingRateConstant_k02.*exp(SensorUnbindingRateConstant_xb2.*F/kT);
        thisPossibleTransition = [thisSensorID thisStateVariableID 0 0 koff F];
        allTransitions(count+1,:) = thisPossibleTransition;
        count = count + 1;
        
    end
    allTransitions(allTransitions(:,1)==0,:) = [];
        
    % [1B] Determine and execute earliest state transition and update time
    % using the Direct Gillespie SSA. [with maxAllowedTimeStep*]
    %       *maxAllowedTimeStep: Impose a max allowed time step. If time to
    %       next state transition exceeds this, advance time (by
    %       maxAllowedTimeStep) and update mechanics without executing a
    %       linkage state transition. This uses the memoryless property of
    %       Markov processes and is implemented to account for
    %       force-dependent rate constants having time dependence.
    % Method: 
    % (1) Determine time to next state transition using inverse
    % transform sampling. tau = 1/sum(rateConstants)*-log(r1) wheere r1 is
    % pseudorandom number on [0,1]
    % (2) If tau <= maxAllowedTimeStep: Determine and execute next state
    % transition based on value of rate constants and using r2, a second
    % pseuodorandom number on [0,1], and advance time by dt=tau. Else:
    % Execute no state transition and advance time by
    % dt=maxAllowedTimeStep.
    
    rateConstants = allTransitions(:,5);
    GillespieActivity = sum(rateConstants); % Activity = sum of transition rates
    r1 = rand();
    dt = -log(r1)/GillespieActivity; % Time to next event. This time is exponentiatlly distributed with rate parameter equal to sum of all rates, and it is computed here via inverse transform sampling.
    
    if dt>maxAllowedTimeStep % Advance time by maxAllowedTimeStep and execute no event
        dt=maxAllowedTimeStep;
        executeTransition = 0;
        numRejectedTransitions = numRejectedTransitions + 1; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt
        
    else % Execute event and advance time by dt from Gillespie SSA
        executeTransition = 1;
        numAcceptedTransitions = numAcceptedTransitions + 1; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt
        % Determine next event.
        probvec = cumsum(rateConstants)/GillespieActivity;
        r2 = rand();
        idNextTransition = find(r2<=probvec,1);
        % Execute next transition
        nextTransition = allTransitions(idNextTransition,:);
        thisSensorID = nextTransition(1);
        thisStateVariableID = nextTransition(2);
        thisStateVariableValue = nextTransition(4);
        if thisStateVariableID==1 % Donor State Transition
            AllSensors_State_Donor(thisSensorID) = thisStateVariableValue;         
        elseif thisStateVariableID==2 % Acceptor State Transition
            AllSensors_State_Acceptor(thisSensorID) = thisStateVariableValue;
        elseif thisStateVariableID==0 % Turnover Transition (D#A#-->D1A1)
            AllSensors_State_Donor(thisSensorID) = 1;
            AllSensors_State_Acceptor(thisSensorID) = 1;
        end          
    end
    
    % Update Time Since Last Turnover for All Sensors
    AllSensors_Lifetime = AllSensors_Lifetime + dt;

    % Update time
    t = t + dt; 
    avgTimeStep = (dt + avgTimeStep*(numRejectedTransitions+numAcceptedTransitions))/(numRejectedTransitions+numAcceptedTransitions+1); % Record simulation performance stats: avg length of time step
    
    % [1C] Optional: Log Linkage Transition 
    if Toggle_OutputLogOfAllTransitions
        if executeTransition && t>timeSpan_LogTransitionData(1) && t<=timeSpan_LogTransitionData(2)
            iDataLog_Transitions = iDataLog_Transitions + 1;
            % Pre-allocate new block if needed
            if iDataLog_Transitions>size(DataLog_Transitions,1)
                DataLog_Transitions = [DataLog_Transitions; PreAllocationBlocks.DataLog_Transitions];
            end
            DataLog_Transitions(iDataLog_Transitions,:) = [t, nextTransition, GillespieActivity, AllSensors_Lifetime(thisSensorID)];
        end
    end
    
    % If Turnover Event Occured
    % Append Lifetime to DataLog_LifetimesAtTurnoverOnly
    % Reset Lifetime to 0
    if thisStateVariableID==0 % Turnover Transition (D#A#-->D1A1)
        DataLog_LifetimesAtTurnoverOnly = [DataLog_LifetimesAtTurnoverOnly; AllSensors_Lifetime(thisSensorID)];
        AllSensors_Lifetime(thisSensorID) = 0;
    end
    
    %% [2] Record Timecourse Data   
    if Toggle_OutputTimeCourseOfEnsemble
        % Log linkage data at regular time intervals
        if (t-t_lastDataLogEnsemble) >= interval_LogEnsembleData
            iDataLog_Ensemble = iDataLog_Ensemble + 1;
            t_lastDataLogEnsemble = t;
            % Pre-allocate new block if needed
            if iDataLog_Ensemble>size(DataLog_Ensemble,1)
                DataLog_Ensemble = [DataLog_Ensemble; PreAllocationBlocks.DataLog_MotorClutch];
            end

            num_D1A1 = sum( AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1 );
            num_D1A0 = sum( AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0 );
            num_D0A1 = sum( AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1 );
            num_D0A0 = sum( AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0 );

            DataLog_Ensemble(iDataLog_Ensemble,:) = [...
                        t,...
                        num_D1A1,...
                        num_D1A0,...
                        num_D0A1,...
                        num_D0A0,...
                        ];
        end
    end
            
end
    
%% Post Process
% Remove unused preallocated storage
DataLog_Transitions(iDataLog_Transitions+1:end,:) = [];
DataLog_Ensemble(iDataLog_Ensemble+1:end,:) = [];
 
% Package up variable names for export
AllVariableNames.Transitions = VariableNames_DataLog_Transitions;
AllVariableNames.Ensemble = VariableNames_DataLog_Ensemble;

num_D1A1 = sum( AllSensors_State_Donor==1 & AllSensors_State_Acceptor==1 );
num_D1A0 = sum( AllSensors_State_Donor==1 & AllSensors_State_Acceptor==0 );
num_D0A1 = sum( AllSensors_State_Donor==0 & AllSensors_State_Acceptor==1 );
num_D0A0 = sum( AllSensors_State_Donor==0 & AllSensors_State_Acceptor==0 );

% Compute Simulation Performance Stats
SimulationPerformanceStats.numAcceptedTransitions = numAcceptedTransitions; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt, so a reaction will occur along with update to forces
SimulationPerformanceStats.numRejectedTransitions = numRejectedTransitions; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt, so no reaction will occur, just update forces
SimulationPerformanceStats.probabilityRejectTransition = numRejectedTransitions/(numAcceptedTransitions+numRejectedTransitions);
SimulationPerformanceStats.avgTimeStep = avgTimeStep; % Record simulation performance stats: avg length of time step

%% XXXXX Support Functions XXXXX
% FP Mechanial Switching, Forward Rate Constant
function out = kms(F,F_half,m,kms0,dxms)
    out = (1./(1+exp(-m*(F-F_half)))) .* kms0.*exp(F.*dxms/kT);
    out = min(out,rateMax);
end

% Protein Unbinding/Unloading
function out = kunload_TwoPathBellModel(F,k01,xb1,k02,xb2,kT)
    out = k01*exp(xb1*F/kT) + k02*exp(xb2*F/kT);
    out = min(out,rateMax);
end

function rate = kBell(F,BellKinetics)
    k0=BellKinetics(1);
    x=BellKinetics(2);
    rate = k0*exp(F*x);
    rate = min(rate,rateMax);
end

end