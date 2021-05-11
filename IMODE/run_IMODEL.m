%% ============ Improved Multi-operator Differential Evolution Algorithm (IMODE) ============
% Should you have any queries, please contact
% Dr. Karam Sallam. Zagazig University
% karam_sallam@zu.edu.eg
% =========================================================================
%% Some part of this code is taken from UMOEA-II
%%==============================================================================
function [outcome,com_time,SR,avgFE,res_det,bestx]= run_IMODEL(run,I_fno,C,D)

Par= Introd_Par(I_fno,D,C);
iter=0;             %% current generation

%% =================== Define a random seed ===============================
%%-----Becase we ran experiments in parallel, we used "*run" to differentiate
%%-----among runs which started at the same time
% RandStream.setGlobalStream (RandStream('mt19937ar','seed',run));
Rand_Seeds=load('input_data/Rand_Seeds.txt');
rand_ind = mod(Par.n*I_fno*30+run-30,length(Rand_Seeds));
run_seed=Rand_Seeds(max(1, rand_ind));
rng(run_seed,'twister');
%% to record seeds for further validation, if needed
% seed_run=stream.Seed;

%% define variables
current_eval=0;             %% current fitness evaluations
PS1=Par.PopSize;            %% define PS1
PS2=100;
InitPop1=PS2;
%% ====================== Initalize x ==================================
x=repmat(Par.xmin,Par.PopSize,1)+repmat((Par.xmax-Par.xmin),Par.PopSize,1).*rand(Par.PopSize,Par.n);

%% calc. fit. and update FES
tic;
fitx= Parametrized_benchmark_func(x',I_fno,C);
current_eval =current_eval+Par.PopSize;
res_det= min(repmat(min(fitx),1,Par.PopSize), fitx); %% used to record the convergence

%% ====================== store the best ==================
[bestold, bes_l]=min(fitx);     bestx= x(bes_l,:);
%% ================== fill in for each  Algorithm ===================================
%% IMODE
EA_1= x(1:PS1,:);    EA_obj1= fitx(1:PS1);   EA_1old = x(randperm(PS1),:);
%% CMA-ES
EA_2= x(PS1+1:size(x,1),:);    EA_obj2= fitx(PS1+1:size(x,1));
%% ================ define CMA-ES parameters ==============================
setting=[]; bnd =[]; fitness = [];
[setting]= init_cma_par(setting,EA_2, Par.n, PS2);

%% ===== prob. of each DE operator
probDE1=1./Par.n_opr .* ones(1,Par.n_opr);
%% ===================== archive data ====================================
arch_rate=2.6;
archive.NP = arch_rate * PS1; % the maximum size of the archive
archive.pop = zeros(0, Par.n); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions
%% ==================== to adapt CR and F =================================
hist_pos=1;
memory_size=20*Par.n;
archive_f= ones(1,memory_size).*0.2;
archive_Cr= ones(1,memory_size).*0.2;
archive_T = ones(1,memory_size).*0.1;
archive_freq = ones(1, memory_size).*0.5;
%%
stop_con=0; avgFE=Par.Max_FES; InitPop=PS1; thrshold=1e-08;

cy=0;
indx=0;
Probs=ones(1,2);
F = normrnd(0.5,0.15,1,PS1);
cr= normrnd(0.5,0.15,1,PS1);

%% main loop
while stop_con==0
    
    iter=iter+1;
    cy=cy+1; % to control CS
    
    %% uncomment if the algorithm uses IMODE and CMA-ES
    %  ================ determine the best phase ===========================
    %     if(cy==ceil(Par.CS+1))
    %
    %         %%calc normalized qualit -- NQual
    %         qual(1) = EA_obj1(1);
    %         qual(2) = EA_obj2(1);
    %         norm_qual = qual./sum(qual);
    %         norm_qual=1-norm_qual; %% to satisfy the bigger is the better
    %
    %         %%Normalized diversity
    %         D(1) = mean(pdist2(EA_1(2:PS1,:),EA_1(1,:)));
    %         D(2) = mean(pdist2(EA_2(2:PS2,:),EA_2(1,:)));
    %         norm_div= D./sum(D);
    %
    %         %%Total Imp
    %         Probs=norm_qual+norm_div;
    %         %%Update Prob_MODE and Prob_CMAES
    %         Probs = max(0.1, min(0.9,Probs./sum(Probs)));
    %
    %         [~,indx]=max(Probs);
    %         if Probs(1)==Probs(2)
    %             indx=0;%% no sharing of information
    %         end
    %
    %
    %     elseif cy==2*ceil(Par.CS)
    %
    %         %% share information
    % %         if indx==1
    % %             list_ind = randperm(PS1);
    % %             list_ind= list_ind(1:(min(PS2,PS1)));
    % %             EA_2(1:size(list_ind,2),:)= EA_1(list_ind,:);
    % %             EA_obj2(1:size(list_ind,2))= EA_obj1(list_ind);
    % %             [setting]= init_cma_par(setting,EA_2, Par.n, PS2);
    % %             setting.sigma= setting.sigma*(1- (current_eval/Par.Max_FES));
    % %         else
    % %             if (min (EA_2(1,:)))> -100 && (max(EA_2(1,:)))<100 %% share best sol. in EA_2 if it is feasible
    % %                 EA_1(PS1,:)= EA_2(1,:);
    % %                 EA_obj1(PS1)= EA_obj2(1);
    % %                 [EA_obj1, ind]=sort(EA_obj1);
    % %                 EA_1=EA_1(ind,:);
    % %             end
    % %
    % %         end
    %         %% reset cy and Probs
    %         cy=1;   Probs=ones(1,2);
    %     end
    Probs = [1 0]; %% comment this line if the algorithm uses IMODE and CMA-ES
    
    %% ======================Applying IMODE ============================
    if (current_eval<Par.Max_FES)
        if rand<=Probs(1)
            
            %% =============================== Linear Reduction of PS1 ===================================================
            UpdPopSize = round((((Par.MinPopSize - InitPop) / Par.Max_FES) * current_eval) + InitPop);
            if PS1 > UpdPopSize
                reduction_ind_num = PS1 - UpdPopSize;
                if PS1 - reduction_ind_num <  Par.MinPopSize
                    reduction_ind_num = PS1 - Par.MinPopSize;
                end
                %% remove the worst ind.
                for r = 1 : reduction_ind_num
                    vv=PS1;
                    EA_1(vv,:)=[];
                    EA_1old(vv,:)=[];
                    EA_obj1(vv)=[];
                    PS1 = PS1 - 1;
                end
                archive.NP = round(arch_rate * PS1);
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            %             end
            
            %% apply IMODE
            [EA_1, EA_1old, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq, current_eval,res_det,F,cr] = ...
                IMODE(EA_1,EA_1old, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,....
                archive_freq, Par.xmin, Par.xmax, Par.n, PS1, current_eval,I_fno,C,res_det,Par.Printing,Par.Max_FES, Par.Gmax, iter,F,cr);
        end
    end
    %% ====================== CMA-ES ======================
    if (current_eval<Par.Max_FES)
        if rand<Probs(2)
            UpdPopSize = round((((Par.MinPopSize1 - InitPop1) / Par.Max_FES) * current_eval) + InitPop1);
            if PS2 > UpdPopSize
                reduction_ind_num = PS2 - UpdPopSize;
                if PS2 - reduction_ind_num <  Par.MinPopSize
                    reduction_ind_num = PS2 - Par.MinPopSize;
                end
                %% remove the worst ind.
                for r = 1 : reduction_ind_num
                    vv=PS2;
                    EA_2(vv,:)=[];
                    %                     EA_1old(vv,:)=[];
                    EA_obj2(vv)=[];
                    PS2 = PS2 - 1;
                end
                %                 archive.NP = round(arch_rate * PS1);
                %                 if size(archive.pop, 1) > archive.NP
                %                     rndpos = randperm(size(archive.pop, 1));
                %                     rndpos = rndpos(1 : archive.NP);
                %                     archive.pop = archive.pop(rndpos, :);
                %                 end
            end
            [ EA_2, EA_obj2, setting,bestold,bestx,bnd,fitness,current_eval,res_det] = ...
                Scout( EA_2, EA_obj2, probSC, setting, iter,bestold,bestx,fitness,bnd,...
                Par.xmin,Par.xmax,Par.n,PS2,current_eval,I_fno,C,res_det,Par.Printing,Par.Max_FES);
            %             disp([min(EA_obj2) Probs(1) Probs(2)]);
        end
    end
    %% ============================ LS2 ====================================
    if current_eval>0.85*Par.Max_FES && current_eval<Par.Max_FES
        if rand<Par.prob_ls
            old_fit_eva=current_eval;
            [bestx,bestold,current_eval,succ] = LS2(bestx,bestold,Par,current_eval,I_fno,C,Par.Max_FES,Par.xmin,Par.xmax);
            if succ==1 %% if LS2 was successful
                EA_1(PS1,:)=bestx';
                EA_obj1(PS1)=bestold;
                [EA_obj1, sort_indx]=sort(EA_obj1);
                EA_1= EA_1(sort_indx,:);
                
                EA_2=repmat(EA_1(1,:), PS2, 1);
                [setting]= init_cma_par(setting,EA_2, Par.n, PS2);
                setting.sigma=1e-05;
                EA_obj2(1:PS2)= EA_obj1(1);
                Par.prob_ls=0.1;
            else
                Par.prob_ls=0.01; %% set p_LS to a small value it  LS was not successful
            end
            %% record best fitness -- set Par.Printing==0 if not
            if Par.Printing==1
                res_det= [res_det repmat(bestold,1,(current_eval-old_fit_eva))];
            end
            
        end
    end
    
    %% ====================== stopping criterion check ====================
    if (current_eval>=Par.Max_FES-4*UpdPopSize)
        stop_con=1;
        avgFE=current_eval;
    end
    if (abs(Par.f_optimal - bestold)<= thrshold)
        stop_con=1;
        bestold=Par.f_optimal;
        avgFE=current_eval;
    end
    
    %% =============================== Print ==============================
    %          fprintf('current_eval\t %d fitness\t %d \n', current_eval, abs(Par.f_optimal-bestold));
    if stop_con
        com_time= toc;%cputime-start_time;
        fprintf('run\t %d, fitness\t %d, avg.FFE\t %d\t %d\n', run, abs(Par.f_optimal-bestold),avgFE,indx(1));
        outcome=abs(Par.f_optimal-bestold);
        if (min (bestx))< -100 || (max(bestx))>100 %% make sure  that the best solution is feasible
            fprintf('in problem: %d, there is  a violation',I_fno);
        end
        SR=(outcome==0);
    end
end
end
