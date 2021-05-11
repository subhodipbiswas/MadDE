%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of LSHADE_cnEpSin which is a new version of LSHADE-EpSin.
%% Please see the following papers:
%% 1. LSHADE_cnEpSin:
%%     Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan, Ensemble Sinusoidal Differential Covariance Matrix Adaptation with Euclidean Neighborhood  for Solving CEC2017 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2017, June, Donostia - San Sebastián, Spain

%% 2. LSHADE-EpSin:
%%    Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan and Robert G. Reynolds: An Ensemble Sinusoidal Parameter Adaptation incorporated with L-SHADE for Solving CEC2014 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2016, Canada, July, 2016

%% About L-SHADE, please see following papers:
%% Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%  J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,” IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945–958, 2009
function LSHADE_cnEpSin(Runs,fhd,CT,problem_size,funcs,max_nfes,optimum)

Rand_Seeds=load('input_data/Rand_Seeds.txt');

Alg_Name=[ 'LSHADE_cnEpSin_(' num2str(CT(1)) num2str(CT(2)) num2str(CT(3)) ')'];

lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];

fprintf('Running %s algorithm on D= %d\n',Alg_Name, problem_size)

for n=0:15
    RecordFEsFactor(n+1) = round(problem_size^((n/5)-3)*max_nfes);
end

progress = numel(RecordFEsFactor);
val_2_reach = 10^(-8);

%%% change freq
freq_inti = 0.5;

pb = 0.4;
ps = .5;

% S.Ndim = problem_size;
% S.Lband = ones(1, S.Ndim)*(-100);
% S.Uband = ones(1, S.Ndim)*(100);
% 
%%%% Count the number of maximum generations before as NP is dynamically
%%%% decreased
G_Max = 0;
if problem_size == 10
    G_Max = 2163;
end
if problem_size == 20
    G_Max = 2745;
end
% if problem_size == 50
%     G_Max = 3022;
% end
% if problem_size == 100
%     G_Max = 3401;
% end

% num_prbs = 30;
% runs = 51;
% allerrorvals = zeros(progress, runs, num_prbs);
% result=zeros(num_prbs,5);

for func_no = funcs
    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func_no, problem_size)
    allerrorvals = zeros(progress, Runs);
    parfor run_id = 1 : Runs
        rand_ind = mod(problem_size*func_no*Runs+run_id-Runs,length(Rand_Seeds));
        run_seed=Rand_Seeds(max(1, rand_ind));
        rng(run_seed,'twister');
        Run_RecordFEsFactor=RecordFEsFactor;
        run_funcvals = [];
        col=1;              %% to print in the first column in all_results.mat
        
        %%  parameter settings for L-SHADE
        EB_best_rate = 0.10;
        arc_rate = 1.4;
        memory_size = 5;
        pop_size = 18 * problem_size;   %18*D
        SEL = round(ps*pop_size);
        
        max_pop_size = pop_size;
        min_pop_size = 4.0;
        
        nfes = 0;
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        
        fitness = feval(fhd,pop',func_no, CT);
        fitness = fitness';
        
        bsf_fit_var = 1e+30;
        bsf_index = 0;
        bsf_solution = zeros(1, problem_size);
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1 : pop_size
            nfes = nfes + 1;
            
            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
                bsf_index = i;
            end
            
            if nfes > max_nfes; break; end
        end
        
        if(nfes>=Run_RecordFEsFactor(1))
            run_funcvals = [run_funcvals; bsf_fit_var];
            Run_RecordFEsFactor(1)=[];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_cr = 0.5 .* ones(memory_size, 1);
        
        memory_freq = freq_inti*ones(memory_size, 1);
        memory_pos = 1;
        
        archive.NP = round(arc_rate * pop_size); % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        
        %% main loop
        gg=0;  %%% generation counter used For Sin
        igen =1;  %%% generation counter used For LS
        
        flag1 = false;
        flag2 = false;
        
        goodF1all = [];
        goodF2all =[];
        badF1all = [];
        badF2all = [];
        goodF1 = [];
        goodF2 = [];
        badF1 = [];
        badF2 = [];
        
        while nfes < max_nfes
            gg=gg+1;
            
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
            
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            mu_freq = memory_freq(mem_rand_index);
            
            %% for generating crossover rate
            cr = normrnd(mu_cr, 0.1);
            term_pos = find(mu_cr == -1);
            cr(term_pos) = 0;
            cr = min(cr, 1);
            cr = max(cr, 0);
            
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            pos = find(sf <= 0);
            
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            
            
            freq = mu_freq + 0.1 * tan(pi*(rand(pop_size, 1) - 0.5));
            pos_f = find(freq <=0);
            while ~ isempty(pos_f)
                freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(length(pos_f), 1) - 0.5));
                pos_f = find(freq <= 0);
            end
            
            sf = min(sf, 1);
            freq = min(freq, 1);
            
            LP = 20;
            flag1 = false;
            flag2 = false;
            if(nfes <= max_nfes/2)
                flag1 = false;
                flag2 = false;
                if (gg <= LP)
                    %% Both have the same probability
                    %% Those generations are the learning period
                    %% Choose one of them randomly
                    p1 = 0.5;
                    p2 = 0.5;
                    c=rand;
                    if(c < p1)
                        sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag1 = true;
                    else
                        sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag2 = true;
                    end
                    
                else
                    %% compute the probability as used in SaDE
                    ns1 = size(goodF1,1);
                    ns1_sum = 0;
                    nf1_sum = 0;
                    %               for hh = 1 : size(goodF1all,2)
                    for hh = gg-LP : gg-1
                        ns1_sum = ns1_sum + goodF1all(1,hh);
                        nf1_sum = nf1_sum + badF1all(1,hh);
                    end
                    sumS1 = (ns1_sum/(ns1_sum + nf1_sum)) + 0.01;
                    
                    
                    ns2 = size(goodF2,1);
                    ns2_sum = 0;
                    nf2_sum = 0;
                    %             for hh = gg-LP : gg-1
                    %               for hh = 1 : size(goodF2all,2)
                    for hh = gg-LP : gg-1
                        ns2_sum = ns2_sum + goodF2all(1,hh);
                        nf2_sum = nf2_sum + badF2all(1,hh);
                    end
                    sumS2 = (ns2_sum/(ns2_sum + nf2_sum)) + 0.01;
                    
                    p1 = sumS1/(sumS1 + sumS2);
                    p2 = sumS2/(sumS2 + sumS1);
                    
                    if(p1 > p2)
                        sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag1 = true;
                        %                   size(goodF1,1)
                    else
                        sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag2 = true;
                        %                   size(goodF2,1)
                    end
                end
            end
            
%             r0 = [1 : pop_size];
            popAll = [pop; archive.pop];

            % Here we select top EDE_best_rate% (10%) to be inculded in EDE
            EDEpNP = max(round(EB_best_rate * pop_size), 2); %% choose at least two best solutions
            EDErandindex = ceil(rand(1, pop_size) .* EDEpNP); %% select from [1, 2, 3, ..., pNP]
            EDErandindex = max(1, EDErandindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            EDEpestind=sorted_index(EDErandindex);
            
            R = Gen_R(pop_size,2);
            R(:,1)=[];
            R=[R EDEpestind];
            fr=fitness(R);
            [~,I] = sort(fr,2);
            R_S=[];
            for i=1:pop_size
                R_S(i,:)=R(i,I(i,:));
            end
            rb=R_S(:,1);
            rm=R_S(:,2);
            rw=R_S(:,3);
            
            vi = pop + sf(:, ones(1, problem_size)) .* (pop(rb, :) - pop + pop(rm, :) - popAll(rw, :));

            vi = boundConstraint(vi, pop, lu);
            
            %       %% Bin Crx
            %       mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            %       rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            %       jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            %       ui = vi; ui(mask) = pop(mask);
            %       %%
            
            %% Bin crossover according to the Eigen coordinate system
            J_= mod(floor(rand(pop_size, 1)*problem_size), problem_size) + 1;
            J = (J_-1)*pop_size + (1:pop_size)';
            crs = rand(pop_size, problem_size) < cr(:, ones(1, problem_size));
            if rand<pb
                %% coordinate ratation
                
                %%%%% Choose neighbourhood region to the best individual
                best = pop(sorted_index(1), :);
                Dis = pdist2(pop,best,'euclidean'); % euclidean distance
                %D2 = sqrt(sum((pop(1,:) - best).^2, 2));
                
                %%%% Sort
                [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                Neighbour_best_pool = pop(idx_ordered(1:SEL), :); %%% including best also so start from 1
                Xsel = Neighbour_best_pool;
                %            sizz = size(Xsel)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %Xsel = pop(sorted_index(1:SEL), :);
                xmean = mean(Xsel);
                % covariance matrix calculation
                C =  1/(SEL-1)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
                C = triu(C) + transpose(triu(C,1)); % enforce symmetry
                [R,D] = eig(C);
                % limit condition of C to 1e20 + 1
                if max(diag(D)) > 1e20*min(diag(D))
                    tmp = max(diag(D))/1e20 - min(diag(D));
                    C = C + tmp*eye(problem_size);
                    [R, D] = eig(C);
                end
                TM = R;
                TM_=R';
                Xr = pop*TM;
                vi = vi*TM;
                %% crossover according to the Eigen coordinate system
                Ur = Xr;
                Ur(J) = vi(J);
                Ur(crs) = vi(crs);
                %%
                ui = Ur*TM_;
                
            else
                
                ui = pop;
                ui(J) = vi(J);
                ui(crs) = vi(crs);
                
            end
            %%%%%%%%
            
            children_fitness = feval(fhd, ui', func_no, CT);
            children_fitness = children_fitness';
            
            
            %%%% To check stagnation
            flag = false;
            bsf_fit_var_old = bsf_fit_var;
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            for i = 1 : pop_size
                nfes = nfes + 1;
                
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i, :);
                    bsf_index = i;
                end
                
                if nfes > max_nfes; break; end
            end
            if(nfes>=Run_RecordFEsFactor(1))
                run_funcvals = [run_funcvals; bsf_fit_var];
                Run_RecordFEsFactor(1)=[];
            end
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            
            dif = abs(fitness - children_fitness);
            
            
            %% I == 1: the parent is better; I == 2: the offspring is better
            I = (fitness > children_fitness);
            goodCR = cr(I == 1);
            goodF = sf(I == 1);
            goodFreq = freq(I == 1);
            dif_val = dif(I == 1);
            
            %% chnage here also
            %% recored bad too
            badF = sf(I == 0);
            
            %% Change Noor
            if flag1 == true
                goodF1 = goodF;
                goodF1all = [goodF1all size(goodF1,1)];
                
                badF1 = badF;
                badF1all = [badF1all size(badF1,1)];
                
                %% Add zero for other one  or add 1 to prevent the case of having NaN
                goodF2all = [goodF2all 1];
                badF2all = [badF2all 1];
                
            end
            if flag2 == true
                goodF2 = goodF;
                goodF2all = [goodF2all size(goodF2,1)];
                
                badF2 = badF;
                badF2all = [badF2all size(badF2,1)];
                
                %% Add zero for other one
                goodF1all = [goodF1all 1];
                badF1all = [badF1all 1];
            end
            %%%%%%
            
            
            %      isempty(popold(I == 1, :))
            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
            [fitness, I] = min([fitness, children_fitness], [], 2);
            
%             run_funcvals = [run_funcvals; fitness];
            
            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);
            
            num_success_params = numel(goodCR);
            
            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;
                
                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                
                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                    memory_cr(memory_pos)  = -1;
                else
                    memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end
                
                %% for updating the memory of freq
                if max(goodFreq) == 0 || memory_freq(memory_pos)  == -1
                    memory_freq(memory_pos)  = -1;
                else
                    memory_freq(memory_pos) = (dif_val' * (goodFreq .^ 2)) / (dif_val' * goodFreq);
                end
                
                memory_pos = memory_pos + 1;
                if memory_pos > memory_size;  memory_pos = 1; end
            end
            
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
            
            if pop_size > plan_pop_size
                reduction_ind_num = pop_size - plan_pop_size;
                if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
                
                pop_size = pop_size - reduction_ind_num;
                SEL = round(ps*pop_size);
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
                
                archive.NP = round(arc_rate * pop_size);
                
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            
        end %%%%%%%%nfes
        
%         bsf_error_val = bsf_fit_var - optimum;
%         if bsf_error_val < val_2_reach
%             bsf_error_val = 0;
%         end

        if(CT(1)==1)
            run_funcvals=run_funcvals-optimum(func_no);
        end
        
        run_funcvals(run_funcvals<val_2_reach)=0;

        fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , run_funcvals(end))
        allerrorvals(:, run_id) = run_funcvals;
        
        %%% Noor: To print files
        %     errorVals= [];
        %     for w = 1 : progress
        %         bestold = run_funcvals(RecordFEsFactor(w) * max_nfes) - optimum;
        %         if abs(bestold)>1e-8
        %             errorVals(w)= abs(bestold);
        % %             col=col+1;
        %         else
        %             bestold=0;
        %             col=col+1;
        %             errorVals(w)= bestold;
        % %             current_eval=max_eval;
        %         end
        %     end
        %     allerrorvals(:, run_id, func) = errorVals;
        
    end %% end 1 run
    
    [~, sorted_index] = sort(allerrorvals(end,:), 'ascend');
    allerrorvals = allerrorvals(:, sorted_index);
    
    fprintf('min_funvals:\t%e\n',min(allerrorvals(end,:)));
    fprintf('median_funvals:\t%e\n',median(allerrorvals(end,:)));
    fprintf('mean_funvals:\t%e\n',mean(allerrorvals(end,:)));
    fprintf('max_funvals:\t%e\n',max(allerrorvals(end,:)));
    fprintf('std_funvals:\t%e\n',std(allerrorvals(end,:)));
    
    file_name=sprintf('Results/%s_%s_%s.txt',Alg_Name,int2str(func_no),int2str(problem_size));
    save(file_name, 'allerrorvals', '-ascii');
    
end %% end 1 function run
end
