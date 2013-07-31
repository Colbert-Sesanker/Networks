function  MALA(Model)

% Start the timer for burn-in samples
tic
% Print variable output to terminal in real time
more off; 
% Note: state variables are called 'species'
% Y is the experimental data we are trying to fit
Y = Model.NoisyData;
[NumOfSpecies, NumTimePts] = size(Y);


% Set Model specifications
% Observed and unobserved species are vectors
% with indexes of state variables
% 'Equations' is a function handle of ODE and
% sensitivity equations
% sampled parameters are the parameters 
% that will be sampled with MCMC
% numSampledParams is the # of parameters
% used in MCMC sampling
numSampledParams                   = Model.numSampledParams;
numTotalParameters                 = Model.numTotalParameters;
SpeciesObserved                    = Model.ObservedSpecies;
SpeciesUnobserved                  = Model.UnobservedSpecies;
totalSpecies                       = length(SpeciesObserved) + ...
                                     length(SpeciesUnobserved);
SDNoiseAdded                       = Model.SDNoiseAdded;
is_AnalyticSens                    = Model.is_AnalyticSens;
if is_AnalyticSens
    Equations                      = Model.Equations_Sens;
    InitialValues                  = Model.InitialValues_Sens;
else
    Equations                      = Model.Equations;
    InitialValues                  = Model.InitialValues;
end
StepSize                           = Model.initialStepSize;
stepSizeRange                      = Model.stepSizeRange; 
Burnin                             = Model.Burnin;
NumOfPosteriorSamples              = Model.NumOfPosteriorSamples;
EquationName                       = Model.EquationName;
totalParameters                    = Model.totalParameters;
TimePoints                         = Model.TimeData;
zero_MetricTensor_derivatives      = Model.zero_MetricTensor_derivatives;
plotProposedTrajectories           = Model.plotProposedTrajectories;

sampled_Param_idxs                 = Model.SampledParameter_idxs;
speciesEstimatesBest               = Model.SpeciesEstimates;
epsilon                            = Model.epsilon;

% Model specific priors:
% Prior struct of priors and derivatives
% Objective priors throw exceptions when parameters reach
% regions of zero probability
Prior                   = Model.Prior;


% Used in calculating prior probabilities of current and proposed parameters
prior                   = Prior.prior;
% Used in computing natural gradient of posterior
prior_derivative        = Prior.prior_derivative;

% Used in computing the metric tensor of posterior
prior_second_derivative = Prior.prior_second_derivative;
% Used in computing the derivative of metric tensor in Laplace
% Beltrami Operator
prior_third_derivative  = Prior.prior_third_derivative;
                          
sampledParameters       = totalParameters(sampled_Param_idxs);

% Set up noise for likelihood function
% Fix noise - CurrentNoise is the variance
% Note: we are cheating here the noise of 
% the likelihood is set to  noise added to the  
% synthetic data
% if SDNoise is a scalar, the same noise is added to each 
% species. Otherwise assume it's a vector of noises for
% each species
if isscalar(SDNoiseAdded),
    CurrentNoise = ones(1, NumOfSpecies) * SDNoiseAdded^2;
end

% Default step size is: 1 /  (number of covariates) ^ (- 1/3)
if strcmp(StepSize, 'default')
    StepSize = numSampledParams^(- 1/3);
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Functions for calculations                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
function [SecondTerm, ThirdTerm] = secondAndThirdTerms(InvG, GDeriv,...
                                                       numSampledParams)  
        for k = 1: numSampledParams
            InvGdG{k}        = InvG * GDeriv{k};
            TraceInvGdG(k)   = trace(InvGdG{k});        
            SecondTerm(:, k) = InvGdG{k}*InvG(:, k);  
        end      
            SecondTerm = sum(SecondTerm, 2)' ;
            ThirdTerm  = InvG * TraceInvGdG' ; 
            ThirdTerm  = ThirdTerm'  ;
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate metric-tensor, gradient, etc. about maximum likelihood        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[  ~    , ...
 trajectories] = ...
          integrate_equations( Equations,     TimePoints,...
                               InitialValues, totalParameters...
                             );

speciesEstimates  = extract_Species_Trajectories(trajectories,...
                                                 NumOfSpecies);
                          
 % calculate 1st and 2nd order sensitivities
if  is_AnalyticSens
    
    [Sensitivities_1,...
     Sensitivities_2] = get_Sensitivities_Analytic(trajectories,... 
                                                   NumOfSpecies,... 
                                                   numSampledParams);
else
    [Sensitivities_1,...
     Sensitivities_2] = ...
                             ...
     get_Sensitivities_FD( Equations,              NumOfSpecies,... 
                           numSampledParams,       totalParameters,...
                           sampled_Param_idxs,     speciesEstimates,...
                           TimePoints,             InitialValues,...
                           epsilon,                zero_MetricTensor_derivatives...                
                         );                           
end % if                                    


gradient_LL = LL_Gradient(...
                            sampledParameters,  Sensitivities_1,... 
                            numSampledParams,   SpeciesObserved,...
                            NumTimePts,         speciesEstimates,... 
                            Y,                  CurrentNoise,...
                            prior_derivative...
                         );
                                                                 
CurrentG = metric_Tensor(...
                          sampledParameters,  Sensitivities_1,... 
                          numSampledParams,   SpeciesObserved,... 
                          CurrentNoise,       prior_second_derivative...                              
                        );

identity           = eye(numSampledParams);    
% In addition to bounding singular values
% add diagonal dust to improve rank                                      
CurrentInvG        = identity / (CurrentG + identity*1e-6);
CurrentFirstTerm   = (CurrentInvG * gradient_LL')' ;

                                    
if zero_MetricTensor_derivatives

   CurrentSecondTerm =  zeros(1, numSampledParams);
   CurrentThirdTerm  =  zeros(1, numSampledParams);
   
else   
   GDeriv   = ...
              metric_Tensor_Deriviatives...
             (... 
               Sensitivities_1,   SpeciesObserved,...
               Sensitivities_2,   numSampledParams,...
               CurrentNoise,      sampledParameters,... 
               prior_third_derivative... 
             );
 
 [CurrentSecondTerm, ...
  CurrentThirdTerm] = secondAndThirdTerms(CurrentInvG, GDeriv,...
                                          numSampledParams);
                  
end % if   
                                   

Current_LL   = calculate_LL( speciesEstimates, Y,... 
                             CurrentNoise,     SpeciesObserved...
                           );                     
                           

currentSampledParams       = sampledParameters;
CurrentSpeciesEstimates    = speciesEstimates;
                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize constants and allocate arrays for Manifold MALA              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up accepted / rejected proposal counters
Accepted  = 0;
Attempted = 0;

% Allocate Histories, LL is for log likelihood
ParaHistory         = zeros(NumOfPosteriorSamples, numTotalParameters);
LL_History          = zeros(NumOfPosteriorSamples, NumOfSpecies);
MetricTensorHistory = cell(1, NumOfPosteriorSamples);
TrajectoryHistory   = cell(1, NumOfPosteriorSamples);

% Set monitor rate for adapting step sizes
MonitorRate         = 10;

% Allocate vector to store acceptance ratios
acceptanceRatios    = zeros(1, Burnin +...
                            NumOfPosteriorSamples);
% Converged to posterior
Converged           = false;
% Set up converged flag
ContinueIterations  = true;
% Initialize iteration number
IterationNum        = 0;
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              MALA Algorithm to sample the parameters:                   %
%          All proposed parameters lack a 'Current' prefix                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ContinueIterations    
    
    IterationNum = IterationNum + 1; 
    
    Attempted    = Attempted    + 1; 
    
    disp(['Iteration:  '...
          num2str(IterationNum)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update parameters       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%          
    
    G     = CurrentG;
    InvG  = CurrentInvG;    
    
    Mean  = currentSampledParams  + (StepSize^2 / 2) *  CurrentFirstTerm  ...
                                  -  StepSize^2      *  CurrentSecondTerm ...
                                  + (StepSize^2 / 2) *  CurrentThirdTerm;
    
    newSampledParams      = Mean +  ...
                            randn(1, numSampledParams) * chol(StepSize^2 * InvG); 
    
    newParams             = update_totalParameters( totalParameters,...
                                                    newSampledParams,...
                                                    sampled_Param_idxs);   
    
    totalMean             = update_totalParameters( totalParameters,...
                                                    Mean,...
                                                    sampled_Param_idxs);  
    
             
    [   ~    , ...
     trajectories]  = ...
               integrate_equations( Equations,     TimePoints,...
                                    InitialValues, newParams...
                                  );                                      
    
    speciesEstimates  = extract_Species_Trajectories(trajectories,... 
                                                     NumOfSpecies);        
           
    % calculate 1st and 2nd order sensitivities
     if  is_AnalyticSens         
         [Sensitivities_1,...
          Sensitivities_2] = get_Sensitivities_Analytic(trajectories,... 
                                                        NumOfSpecies,... 
                                                        numSampledParams);
     else
         [Sensitivities_1,...
          Sensitivities_2] = ...
                             ...
          get_Sensitivities_FD( Equations,              NumOfSpecies,... 
                                numSampledParams,       newParams,...
                                sampled_Param_idxs,     speciesEstimates,...
                                TimePoints,             InitialValues,...
                                epsilon,                zero_MetricTensor_derivatives...                
                              );                           
     end % if   
                                   
         
    % plot trajectories as proposed 
    if  plotProposedTrajectories
        for i = 1:NumOfSpecies
            figure(i);
            hold on;
            plot(speciesEstimates(i, :));
        end % for
    end % if
    
    % Calculate probability of proposed parameters given current parameters
    % The first term is the log of the normalization constant
    ProbNewGivenOld   = - log(prod(diag(chol(InvG * StepSize^2)))) - ...
                          0.5*(Mean - newSampledParams) * (G / StepSize^2) * (Mean - newSampledParams)'... 
                          -   (numSampledParams / 2)*log(2*pi);   
    try
        gradient_LL = LL_Gradient(...
                                   newSampledParams,   Sensitivities_1,... 
                                   numSampledParams,   SpeciesObserved,... 
                                   NumTimePts,         speciesEstimates,... 
                                   Y,                  CurrentNoise,...
                                   prior_derivative...
                                 );
    catch       
        IterationNum = IterationNum - 1;    
        Attempted    = Attempted    - 1;  
        % redo current iteration step 
        continue
    end % try                                                                                                     
                                  
   
    G           = metric_Tensor(...
                                 newSampledParams,   Sensitivities_1,... 
                                 numSampledParams,   SpeciesObserved,... 
                                 CurrentNoise,       prior_second_derivative...                              
                               );
                               
    identity    = eye(numSampledParams); 
    % In addition to bounding singular values
    % add diagonal dust to improve rank                                       
    InvG        = identity / (G + identity*1e-6);
    FirstTerm   = (CurrentInvG * gradient_LL')';  
    
    if zero_MetricTensor_derivatives

     SecondTerm =  zeros(1, numSampledParams);
     ThirdTerm  =  zeros(1, numSampledParams);
   
    else     
        GDeriv = ...
                 metric_Tensor_Deriviatives...
                  (... 
                    Sensitivities_1,   SpeciesObserved,  ...
                    Sensitivities_2,   numSampledParams, ...
                    CurrentNoise,      newSampledParams, ...
                    prior_third_derivative...
                  );  
                                
     [SecondTerm, ...
      ThirdTerm] = secondAndThirdTerms(InvG, GDeriv,...
                                       numSampledParams);
    end % if  
        
        
    Mean  = newSampledParams     + (StepSize^2 / 2) * FirstTerm  ...
                                 -  StepSize^2      * SecondTerm ...
                                 + (StepSize^2 / 2) * ThirdTerm;
                          
    totalMean             = update_totalParameters( totalParameters,...
                                                    Mean,...
                                                    sampled_Param_idxs);                      
     
    % Calculate probability of current parameters given proposed parameters   
    % The first term is the log of the normalization constant 
    % log(prod(diag(chol... is a fast calculation of log determinant                
    ProbOldGivenNew  = - log(prod(diag(chol(InvG * StepSize^2)))) - ...                         
                         0.5*(Mean - currentSampledParams) * (G / StepSize^2) * (Mean - currentSampledParams)'... 
                         -   (numSampledParams / 2)*log(2*pi); 
                         
    Proposed_LL      = calculate_LL( speciesEstimates, Y,... 
                                     CurrentNoise,     SpeciesObserved...
                                   );
    
    % Calculate the log prior for proposed parameter value   
    ProposedLogPrior = zeros(1, numSampledParams);
    for p = 1: numSampledParams    
        ProposedLogPrior(p) = prior(sampled_Param_idxs(p),...
                                    newSampledParams(p));        
    end
    ProposedLogPrior = sum(ProposedLogPrior);                                                                 
                          
    
    % Calculate the log prior for current hyperparameter value
    CurrentLogPrior = zeros(1, numSampledParams);
    for p = 1: numSampledParams    
        CurrentLogPrior(p) = prior(sampled_Param_idxs(p),...
                                   currentSampledParams(p));        
    end
    CurrentLogPrior = sum(CurrentLogPrior);  
        
   % Accept according to ratio of log probabilities
    ratio =... 
                 Proposed_LL  +  ProposedLogPrior  +  ProbOldGivenNew -... 
                 Current_LL   -  CurrentLogPrior   -  ProbNewGivenOld;           
      
        
    if ratio > 0 || log(rand) < min(0, ratio)
        % Accept proposal
        % Update variables        
        Accepted                   = Accepted + 1;
        
        currentParams              = newParams;
        currentSampledParams       = newSampledParams;
        Current_LL                 = Proposed_LL;            
        CurrentG                   = G;
        CurrentInvG                = InvG;            
        CurrentFirstTerm           = FirstTerm;
        CurrentSecondTerm          = SecondTerm;
        CurrentThirdTerm           = ThirdTerm;  
        CurrentSpeciesEstimates    = speciesEstimates;
                   
    end
    
    % Keep track of acceptance ratios at each iteration number
    acceptanceRatio                = 100*Accepted / Attempted;
    acceptanceRatios(IterationNum) = acceptanceRatio; 
    
    disp(['acceptance probability:  ' num2str(ratio) ]);    
    disp(['acceptance ratio:  '       num2str(100*Accepted / Attempted) ]);
    disp(['Parameters:  '             num2str(newSampledParams) ]);
    
    figStart = NumOfSpecies + 1;
    figEnd   = figStart + numSampledParams;    
    plotTraces(currentSampledParams, ...
               figStart: figEnd,     ...
               IterationNum);
    
    cp = currentSampledParams;   
    paramPairs = {... 
                   [cp(1) cp(2)],... 
                   [cp(1) cp(3)],...
                   [cp(2) cp(3)] ...
                 };    
    figStart = figEnd + 1;
    figEnd   = figStart + length(paramPairs);               
    plot_2D_Slice( paramPairs, ...
                   figStart: figEnd...
                 );
    
    paramTriples = {... 
                     [cp(1) cp(2) cp(3)],...                      
                   }; 
                     
    figStart = figEnd + 1;
    figEnd   = figStart + length(paramTriples);
    plot_3D_Slice( paramTriples, ...
                   figStart: figEnd...
                 );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save parameters, LL, metric tensors and trajectories %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Converged
        
        posteriorSampleNum = IterationNum -... 
                             Burnin;
        
        ParaHistory(posteriorSampleNum, :)         = currentParams;
        LL_History(posteriorSampleNum, :)          = Current_LL;        
        MetricTensorHistory{posteriorSampleNum}    = CurrentG;
        TrajectoryHistory{posteriorSampleNum}      = CurrentSpeciesEstimates;        
        
        if IterationNum == Burnin + NumOfPosteriorSamples
            % N Posterior samples have been collected so stop
            ContinueIterations = false;
        end     
        

    else        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust step size based on acceptance ratio  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        if mod(IterationNum, MonitorRate) == 0            
            minAccept = stepSizeRange(1); 
            maxAccept = stepSizeRange(2);
            % amount to increase/decrease step size 
            % to alter acceptance ratio 
            dec       = 0.05; % Note: this should be based on acceptanceRatio - minAccept
            assert(StepSize > 0);            
            
            if     acceptanceRatio < minAccept
                       StepSize = StepSize  - dec;
            elseif acceptanceRatio > maxAccept  
                % Steps are too small, so increase them a bit
                       StepSize = StepSize  + dec;                         
            end % if   
            
            if StepSize < dec 
               StepSize = dec;
            end
            
            disp('%%%%');
            disp(['step size: ' num2str(StepSize) ]);
            disp('%%%%');
                                 
        end % if            
        
        % Change converged tab if burn-in complete
            
        if IterationNum == Burnin          
            Converged  = true;    
            % End burn-in-samples timer: Get time it took to sample Burn-In
            BurnInTime = toc;
            % Begin posterior samples timer: Restart timer to get time for collecting posterior samples
            tic;            
        end % if        
        
    end % if
     
end % while

% End posterior samples timer: Time to collect Posterior Samples
PosteriorTime = toc;

% Save posterior
fileName = [ 'mMALA'...
             '_'... 
             EquationName...              
           ];

save(  ['./Results/' fileName],... 
       'ParaHistory',... 
       'MetricTensorHistory',... 
       'LL_History',... 
       'BurnInTime',... 
       'PosteriorTime',... 
       'Y',... 
       'TimePoints',...
       'acceptanceRatios',...
       'TrajectoryHistory',...
       'StepSize'...
    );


end % main function


