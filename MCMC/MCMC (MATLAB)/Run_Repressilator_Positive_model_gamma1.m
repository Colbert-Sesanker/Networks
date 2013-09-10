% Run Repressilator with positive Feedback

% Add all dierectories in MCMC methods to path
addpath(genpath('./'))
% Close all Figures
close all;

Model.burnin                = 50;
Model.numPosteriorSamples   = 350;


% name for saving results
Model.equationName                   = 'Repressilator_Positive_gamma1'; 
% function handle of model equations 
Model.equations                      = @RepressilatorPositive;
% function handle of model for automatic differentiation
Model.equations_AD                   = @RepressilatorPositive_AD;
% indexes of observed species in state vector
Model.observedStates                 = [1 2 3];
Model.unobservedStates               = [];
Model.totalStates                    = 1:3;
% indexes of sampled parameters 
% All other parameters are held fixed
sampledParam_idxs                    = [1 5 9];
Model.numSampledParams               = length(sampledParam_idxs);
Model.sampledParam_idxs              = sampledParam_idxs;
% Noise
Model.addedNoise_SD                  = 0.5;
% The initial step size for parameter updates
Model.initialStepSize                = .75;
% Step Size for standard Metropolis Hastings
Model.mhStepSize                     = 0.7;
% The step size is adjusted online until acceptance ratio
% is in the range given by 'stepSizeRange'
Model.stepSizeRange                  = [70 80];
% Adjust Step-Size after stepSizeMonitorRate iterations
Model.stepSizeMonitorRate            = 20;
% epsilon is for finite differences
Model.epsilon                        = 5e-1; 
Model.zeroMetricTensorDerivatives    = true;
% If true plots trajectories for all proposed parameters
Model.plotProposedTrajectories       = true;
% Use basic MALA algorithm without manifold information
Model.isMala                         = false;
% Number of steps to recalculate metric tensor (can set to 'randomWalk')
Model.tensorMonitorRate              = 1;
% Preconditioning matrix that don't use 
Model.preConditionMatrix             = eye(Model.numSampledParams);
%                                        inv(1e4*[4.0493   -5.3057   -0.8575
%                                        -5.3057    8.3806    1.5185
%                                        -0.8575    1.5185    0.2977]);

% For RMHMC
Model.numLeapFrogSteps               = 3;
Model.stepSize_RMHMC                 = 2 / Model.numLeapFrogSteps;
Model.maxFixedPointStepsMomentum     = 5;
Model.maxFixedPointStepsPosition     = 3;
% Use basic HMC without manifold information
Model.isHMC                          = false;

% Choose sensitvity type
sensitivityMethods                   = getSensitivityMethods();
% 1 = symbolic, 2= finite diff, 3= automatic diff
Model.sensitivityMethod              = sensitivityMethods{3};


% Model specific priors (function handles))
% The Prior Struct holds all information regarding priors

% Used in calculating prior probabilities of current and proposed parameters
Prior.prior                   = @gammaPrior;
% Used in computing natural gradient of posterior
Prior.priorDerivative         = @gammaPriorDeriv;
% Used in computing the metric tensor of posterior
shape                         = 1;
Prior.priorSecondDerivative   = @(numParams, sampledParams) ...
                                  gammaPriorSecondDeriv(numParams, sampledParams, shape);
                              
% Used in computing the derivative of metric tensor for the Laplace
% Beltrami Operator
Prior.priorThirdDerivative    = @repressilatorPriorThirdDerivative;

Model.Prior                   = Prior;

% Model specific integration settings
numTimePts    = 200;
startTime     = 0;
endTime       = 200;
step          = (endTime - startTime) / numTimePts;
timePoints    = startTime: step :endTime;

A_initial =     .48199;
B_initial =    5.11385; 
C_initial =  105.77422;

initialValues = [ A_initial...
                  B_initial...
                  C_initial ];         

% All Repressilator parameters 
k1 = 2.35804;   k2 = 4.42269;   k3 = 4.80922;   k4 = 5;
n1 = 5.03005;   n2 = 5.73448;   n3 = 6.05167;   n4 = 7;
a1 = 5.73702;   a2 = 6.92108;   a3 = 7.46407;   a4 = .7;
b1 =  .3284;    b2 =  .4967;    b3 =  .4518;
y1 =  .908;     y2 =  .8093;    y3 = 1.1444;

params =     [...
               k1 k2 k3 k4 ...               
               n1 n2 n3 n4 ...
               a1 a2 a3 a4 ...
               b1 b2 b3    ...
               y1 y2 y3    ...              
             ];

% Cell array of all parameter names         
paramMap = {...
               'k1' 'k2' 'k3' 'k4' ...               
               'n1' 'n2' 'n3' 'n4' ...
               'a1' 'a2' 'a3' 'a4' ...
               'b1' 'b2' 'b3'      ...
               'y1' 'y2' 'y3'      ...              
           };
           
stateMap = {'A' 'B' 'C'};
    
Model.totalParams        = params;
Model.numTotalParams     = length(params);
Model.paramMap           = paramMap;
Model.stateMap           = stateMap;
Model.initialValues      = initialValues;

% Integrate model equations
[timeData, speciesEstimates] = ode45( @RepressilatorPositive,...
                                      timePoints,...
                                      initialValues,...      
                                      odeset('RelTol', 1e-6),...
                                      params);

Model.timeData          = timeData' ; % note the transposes here
speciesEstimates        = speciesEstimates' ;
Model.speciesEstimates  = speciesEstimates;

[numStates, numTimePts] = size(speciesEstimates);

% Add Noise to trajectories
addedNoise_SD    = Model.addedNoise_SD;
Model.noisyData  = speciesEstimates + ...
                   randn(numStates, numTimePts) .* addedNoise_SD;

%%%%%%%%%%%%%%%%%%%%%%%%%        
% Call sampling routines %
%%%%%%%%%%%%%%%%%%%%%%%%%


% MH(Model);
% MH_oneParamAt_a_Time(Model);
  MALA(Model);
% RMHMC(Model);

%%%%%%%%%%%%%%%%%%%%%%%%%        
% For ensemble analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%

% 95%, 50% and 5% confidence intervals on trajectories
Model.trajQuantiles                    = [.975 .5 .025];
Model.paramQuantiles.params            = Model.sampledParam_idxs;
Model.trajectoryQuantiles.statesToPlot = Model.totalStates;
Model.posteriorParamsToPlot            = Model.sampledParam_idxs;
ensembleData                           = load(['mMALA_' Model.equationName]);
ensembleAnalysis(Model, ensembleData);















