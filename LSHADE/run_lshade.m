function run_lshade(Runs,fhd,C,problem_size,funcs,max_nfes,pop_size,optimum)

Rand_Seeds=load('input_data/Rand_Seeds.txt');

Alg_Name=[ 'LSHADE_(' num2str(C(1)) num2str(C(2)) num2str(C(3)) ')'];

F =0.50*ones(pop_size,problem_size);
Cr=0.90*ones(pop_size,problem_size);

lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];

fprintf('Running %s algorithm on D= %d\n',Alg_Name, problem_size)

for n=0:15
    RecordFEsFactor(n+1) = round(problem_size^((n/5)-3)*max_nfes);
end

progress = numel(RecordFEsFactor);
val_2_reach = 10^(-8);

for func_no = funcs
    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func_no, problem_size)
    allerrorvals = zeros(progress, Runs);
%     you can use parfor if you have MATLAB Parallel Computing Toolbox
    parfor run_id = 1 : Runs
        rand_ind = mod(problem_size*func_no*Runs+run_id-Runs,length(Rand_Seeds));
        run_seed=Rand_Seeds(max(1, rand_ind));
        rng(run_seed,'twister');
        Run_RecordFEsFactor=RecordFEsFactor;
        run_funcvals = [];
        
        %%  parameter settings for L-SHADE
        p_best_rate = 0.11;
        arc_rate = 1.4;
        memory_size = 5;
        pop_size = 18 * problem_size;
        
        max_pop_size = pop_size;
        min_pop_size = 4.0;
        
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        
        fitness = feval(fhd, pop', func_no, C);
        fitness = fitness';
        
        nfes = 0;
        bsf_fit_var = 1e+30;
        bsf_solution = zeros(1, problem_size);
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1 : pop_size
            nfes = nfes + 1;
            
            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
            end
            
            if nfes > max_nfes
                break;
            end
        end
        
        if(nfes>=Run_RecordFEsFactor(1))
            run_funcvals = [run_funcvals; bsf_fit_var];
            Run_RecordFEsFactor(1)=[];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_cr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;

        archive.NP = round(arc_rate * pop_size); % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        %% main loop
        while nfes < max_nfes
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
            
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            
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
            
            sf = min(sf, 1);
            
            r0 = [1 : pop_size];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
            
            pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
            
            vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop, lu);
            
            mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            
            children_fitness = feval(fhd, ui', func_no ,C);
            children_fitness = children_fitness';
            
            for i = 1 : pop_size
                nfes = nfes + 1;
                
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                end
                
                if nfes > max_nfes
                    break;
                end
            end
            
            if(nfes>=Run_RecordFEsFactor(1))
                run_funcvals = [run_funcvals; bsf_fit_var];
                Run_RecordFEsFactor(1)=[];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            
            dif = abs(fitness - children_fitness);
            
            %% I == 1: the offspring is better; I == 0: the parent is better
            I = (fitness > children_fitness);
            goodCR = cr(I == 1);
            goodF = sf(I == 1);
            dif_val = dif(I == 1);
               
            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));            
            [fitness, I] = min([fitness, children_fitness], [], 2);
            
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
                
                memory_pos = memory_pos + 1;
                
                if memory_pos > memory_size
                    memory_pos = 1;
                end
            end
            
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
            
            if pop_size > plan_pop_size
                reduction_ind_num = pop_size - plan_pop_size;
                if pop_size - reduction_ind_num <  min_pop_size
                    reduction_ind_num = pop_size - min_pop_size;
                end
                
                pop_size = pop_size - reduction_ind_num;
                
                for r = 1 : reduction_ind_num
                    [valBest, indBest] = sort(fitness, 'ascend');
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
        end
        
        if(C(1)==1)
            run_funcvals=run_funcvals-optimum(func_no);
        end
        
        run_funcvals(run_funcvals<val_2_reach)=0;
        
        fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , run_funcvals(end))
        allerrorvals(:, run_id) = run_funcvals;
        
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
