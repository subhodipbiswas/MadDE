%% ============ Improved Multi-operator Differential Evolution Algorithm (IMODE) ============
% Should you have any queries, please contact
% Dr. Karam Sallam. Zagazig University
% karam_sallam@zu.edu.eg
% =========================================================================
mex cec21_basic_func.cpp -DWINDOWS
mex cec21_bias_func.cpp -DWINDOWS
mex cec21_bias_rot_func.cpp -DWINDOWS
mex cec21_bias_shift_func.cpp -DWINDOWS
mex cec21_bias_shift_rot_func.cpp -DWINDOWS
mex cec21_rot_func.cpp -DWINDOWS
mex cec21_shift_func.cpp -DWINDOWS
mex cec21_shift_rot_func.cpp -DWINDOWS

clc;
clear all;
format short e;
%%  introductory Definitions
num_prbs=10;                        %% number of test problems
max_runs=30;                        %% number of runs
outcome=zeros(max_runs,1);          %% to save the solutions of each run
com_time=zeros(max_runs,1);         %% Computational time
SR=zeros(max_runs,1);               %% How many times the optimal solution is obtained
Avg_FES=zeros(max_runs,1);          %% average fitness evaluations
Final_results=zeros(num_prbs,8);    %% to save the final results

%% Benchmark info

% C is the Parametrized Selection Vector, where:
% C1 is Shift indicator (0=0 ,1=F*)
% C2 Translation indicator (0=0, 1=oi)
% C3 Rotation indicator (0=I 1=M)

C = [
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1;
    ];

%% run on more than one processor
myCluster = parcluster('local');
myCluster.NumWorkers = 10;  % define how many processors to use

%% ========================= main loop ====================================
for D= [10 20]
    
    for i=1:size(C,1)
        Alg_Name=[ 'IMODE_(' num2str(C(i,1)) num2str(C(i,2)) num2str(C(i,3)) ')'];
        
        for I_fno= 1:10
            Par= Introd_Par(I_fno,D,C(i,:)); %% set of parameters
            sol=zeros(max_runs,Par.n); %% the best solution vector of each run
            vv=[];
            
            parfor run=1:max_runs
                [outcome(run),com_time(run),SR(run),Avg_FES(run),res,sol(run,:)]=run_IMODEL(run,I_fno,C(i,:),D);

                %% to print the convergence of ech run % set 0 if not
                if Par.Printing==1
                    res= res-repmat(Par.f_optimal,1,size(res,2));
                    res(res<=1e-08)=0; ss=size(res,2);
                    endv=res(ss);
                    if size(res,2)<Par.Max_FES
                        res(size(res,2):Par.Max_FES)=endv;
                    end
                    vv(run,:)= res(1:Par.Max_FES);
                end
            end
            %% save the results in a text
            % Final_results(I_fno,:)= [min(outcome),max(outcome),median(outcome), mean(outcome),std(outcome),mean(com_time),mean(SR),mean(Avg_FES)];
            % disp(Final_results);
            % save('results.txt', 'Final_results', '-ascii');

            %% fitness values at different levels of the optimization process
            %%% required by the competition
            if Par.Printing==1
                for k=1:16
                    lim(k)=Par.n^(((k-1)/5)-3).*Par.Max_FES;
                end
                lim= ceil(lim);
                % lim= [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].*Par.Max_FES;
                res_to_print= vv(:,lim);
                res_to_print=res_to_print';

                [~, sorted_index] = sort(res_to_print(end,:), 'ascend');
                res_to_print = res_to_print(:, sorted_index);

                fprintf('min_funvals:\t%e\n',min(res_to_print(end,:)));
                fprintf('median_funvals:\t%e\n',median(res_to_print(end,:)));
                fprintf('mean_funvals:\t%e\n',mean(res_to_print(end,:)));
                fprintf('max_funvals:\t%e\n',max(res_to_print(end,:)));
                fprintf('std_funvals:\t%e\n',std(res_to_print(end,:)));
    
                file_name=sprintf('Results/%s_%s_%s.txt',Alg_Name,int2str(I_fno),int2str(D));
                save(file_name, 'res_to_print', '-ascii');

            end
        end
    end
end
