%% ************************************************************************
% About: Surrogate-based Bayesian Hyperparameter Optimizer (SUBHO).
% Author: Subhodip Biswas, Debanjan Saha, Shuvodeep De, Adam D Cobb, Swagatam Das and Brian Jalaian
% Cite us: S. Biswas, D. Saha, S. De, A. D. Cobb, S. Das and B. A. Jalaian, 
% "Improving Differential Evolution through Bayesian Hyperparameter Optimization," 
% 2021 IEEE Congress on Evolutionary Computation (CEC), Kraków, Poland, 2021.
%
% For any queries please feel free to contact us
% Subhodip Biswas : sub17was(at)gmail.com
% Debanjan Saha   : debanjansh(at)gmail.com
mex cec21_basic_func.cpp -DWINDOWS
mex cec21_bias_func.cpp -DWINDOWS
mex cec21_bias_rot_func.cpp -DWINDOWS
mex cec21_bias_shift_func.cpp -DWINDOWS
mex cec21_bias_shift_rot_func.cpp -DWINDOWS
mex cec21_rot_func.cpp -DWINDOWS
mex cec21_shift_func.cpp -DWINDOWS
mex cec21_shift_rot_func.cpp -DWINDOWS
mkdir('Hyperparameters')
mkdir('Results')
mkdir('Solutions')
pause(5);
clear;
clc;
%% Hyper parameter space
% Define the optimization variables corresponding the hyperparameter space
% In this case, we are using Mad DE hyper-parameters:
% probability of doing q-best crossover:  q_cr_rate [0.01, 0.05];
% selecting the top p% of solutions:   p_best_rate [0.05, 0.25];
% archive rate multiplier:     arc_rate [1.0, 3.0];
% memory_size multiplier:      mem_mult [1, 10]; 
% pop_size multiplier:         pop_mult [1.0, 5.0]; 
% initial scaling factor       sf_init  [0.1, 0.9];    
% initial crossover rate       cr_init  [0.1, 0.9];
num_hyper = 7;   % numbr of parameter to tune

%% Design space
% Define the design variables for bayesopt to operate on. Here we use
% integers so that we may control the MadDE hyperparamters with 2 decimal
% places.
for i=1:num_hyper
    switch i
        case 1
            varname = 'q_cr_rate';
            var(i)= optimizableVariable(varname, [1, 5], 'Type', 'integer');
        case 2
            varname = 'p_best_rate';
            var(i)= optimizableVariable(varname, [5, 25], 'Type', 'integer');
        case 3
            varname = 'arc_rate';
            var(i)= optimizableVariable(varname, [100, 300], 'Type', 'integer');
        case 4
            varname = 'mem_mult';
            var(i)= optimizableVariable(varname, [100, 1000], 'Type', 'integer');
        case 5
            varname = 'pop_mult';
            var(i)= optimizableVariable(varname, [100, 500], 'Type', 'integer');
        case 6
            varname = 'sf_init';
            var(i)= optimizableVariable(varname, [10, 90], 'Type', 'integer');
        case 7
            varname = 'cr_init';
            var(i)= optimizableVariable(varname, [10, 90], 'Type', 'integer');
        otherwise
            disp('Error..')
    end        
end

%% Define the objective function handle
fun = @(var) tune_hyperparameters(var, num_hyper);

%% Tune MadDE using SUrrogate-based Bayesian Hyperparameter Optimizer.
% Run the bayesopt for 'budget' FEs and get the solution
rng(1);   % you can use any seed other than 1
budget = 128;  % Define the maximum number of function evaluations (FEs)
warm_start = true; % use this if you already have some stored hyperparamter values
fprintf('\n\n Starting SUBHO at %s.\n\n', string(datetime));

if warm_start
    prev_parms = readHPTs('Hyperparameters');
    [numHPTs, ~] = size(prev_parms);
    if numHPTs > 0
        rec_parms = table(prev_parms(:,1),prev_parms(:,2),prev_parms(:,3),...
            prev_parms(:,4),prev_parms(:,5),prev_parms(:,6),prev_parms(:,7));

        score_1 = prev_parms(:,8);
        score_2 = prev_parms(:,9);
        objective = score_2 + 100 - score_1;
        % We assume that MadDE has already been run using manually-found
        % hyperparameters. So we run bayesopt to find new hyperparameters.
        results = bayesopt(fun, var, 'Verbose', 1,'UseParallel', false,...
            'AcquisitionFunctionName', 'expected-improvement',...
            'IsObjectiveDeterministic', false,...
            'MaxObjectiveEvaluations', numHPTs + budget,...
            'InitialX', rec_parms,...
            'InitialObjective', objective,...
            'SaveVariableName','BayesIterations',...
            'OutputFcn',[],'PlotFcn',[]);
    else
        warm_start = false;
    end
end

if ~warm_start
    % Run MadDE from scratch to with the manually found hyperparameters
    % MadDE(q_cr_rate, p_best_rate, arc_rate, mem_mult, pop_mult, sf_init, cr_init, 'MadDEv0.0');    
    fprintf('\n Generating results of MadDE with manually-found hyperparameters.\n')
    MadDE(1, 10, 140, 500, 250, 50, 50, 'MadDEv0.0');
    fprintf('\n Proposing automated hyperparameters using BayesOpt...\n\n')
    % use bayesopt to find better hyperparameters
    results = bayesopt(fun, var,...
        'Verbose', 1,...
        'UseParallel', false,...
        'AcquisitionFunctionName', 'expected-improvement',...
        'IsObjectiveDeterministic', true,...
        'MaxObjectiveEvaluations', budget,...
        'OutputFcn',[],'PlotFcn',[]);
end

%% Get the results
clc;
fprintf('\n The best hyperparameters found by SUBHO are as follows...\n')
bayesian_hpts = table2array(results.XAtMinEstimatedObjective);
disp(bayesian_hpts); % Print the result of the run
% Call MadDE using using values stored in 'bayesian_hpts'
q_cr_rate = bayesian_hpts(1);
p_best_rate = bayesian_hpts(2);
arc_rate = bayesian_hpts(3);
mem_mult = bayesian_hpts(4);
pop_mult = bayesian_hpts(5);
sf_init = bayesian_hpts(6);
cr_init = bayesian_hpts(7);
fprintf('\n Generating results of MadDE with hyperparameters found by SUBHO...\n')
MadDE(q_cr_rate, p_best_rate, arc_rate, mem_mult, pop_mult, sf_init,...
    cr_init, 'MadDEvHPT')
% Compute the scores of MadDE with the two parameteric configuration
[score_1, score_2] = getScoreopt();
fprintf('\n\nScores of MadDE with \n');
fprintf('\t      manually-found hyperparameters: %2.2f\n', score_2);
fprintf('\t automatically-found hyperparameters: %2.2f\n', score_1);