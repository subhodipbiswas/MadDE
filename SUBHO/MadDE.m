function MadDE(q_cr, p_best, arc_rate, mem_mult, pop_mult, sf_init, cr_init, alg_name)

%%  parameter settings for MadDEvHPT
q_cr_rate = q_cr/100; 
p_best_rate = p_best/100; 
arc_rate = arc_rate/100; 
mem_mult = mem_mult/100; 
pop_mult = pop_mult/100; 
sf_init= sf_init/100;
cr_init = cr_init/100;

format long;
num_runs=30;
fhd=@Parametrized_benchmark_func;
p = parpool();    % Starting the parallel pool
Rand_Seeds=load('input_data/Rand_Seeds.txt');
optimum= [100, 1100, 700, 1900, 1700, 1600, 2100, 2200, 2400, 2500];

%% C is the Parametrized Selection Vector, where:
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

for problem_size = [10 20]
    switch problem_size
        case 10
            max_nfes = 200000;
        case 20
            max_nfes = 1000000;
        otherwise
            disp('Error..')
    end
    
    lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];

    for n=0:15
        RecordFEsFactor(n+1) = round(problem_size^((n/5)-3)*max_nfes);
    end

    progress = numel(RecordFEsFactor);  
    val_2_reach = 10^(-8);
    
    for cc=1:size(C,1)
        Alg_Name=[ alg_name '_(' num2str(C(cc,1)) num2str(C(cc,2)) num2str(C(cc,3)) ')'];

        for func_no = 1:10
            allerrorvals = zeros(progress, num_runs);
            allsolutions = zeros(problem_size, num_runs);

            %%     you can use parfor if you have MATLAB Parallel Computing Toolbox
            parfor run_id = 1 : num_runs
                  rand_ind = mod(problem_size*func_no*num_runs+run_id-num_runs, length(Rand_Seeds));
                  run_seed=Rand_Seeds(max(1, rand_ind));
                  rng(run_seed);
                  Run_RecordFEsFactor=RecordFEsFactor;
                  run_funcvals = [];

                  memory_size = round(mem_mult * problem_size);
                  pop_size = round(pop_mult * problem_size ^ 2);
                  
                  max_pop_size = pop_size;
                  min_pop_size = 4.0;

                  %% Initialize the main population
                  % https://www.mathworks.com/help/stats/sobolset.html
                  % p = sobolset(problem_size, 'Skip', 1e4, 'Leap', 1e3);
                  % p = scramble(p, 'MatousekAffineOwen');
                  % rand0 = net(p, pop_size);
                  rand0 = rand(pop_size, problem_size);
                  popold = repmat(lu(1, :), pop_size, 1) + rand0 .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
                  pop = popold; % the old population becomes the current population

                  nfes = 0;                               %% current evaluation
                  bsf_fit_var = 1e+30;                    %% best fitness
                  bsf_solution = zeros(1, problem_size);
                  
                  %%  Calculate functional evaluation on CEC21 Benchmark Suite
                  fitness = feval(fhd, pop', func_no, C(cc,:));
                  fitness = fitness';

                  %%%%%%%%%% record and store fitness  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

                  %% ===== prob. of each DE operatorr =
                  num_de = 3;
                  count_S = zeros(1, num_de);
                  probDE = 1./num_de .* ones(1, num_de);

                  memory_sf = sf_init .* ones(memory_size, 1);  % 0.5 .* ones(memory_size, 1);
                  memory_cr = cr_init .* ones(memory_size, 1);  % 0.5 .* ones(memory_size, 1);
                  memory_pos = 1;

                  archive.NP = round(arc_rate * pop_size); % the maximum size of the archive
                  archive.pop = zeros(0, problem_size); % the solutions stored in te archive
                  archive.funvalues = zeros(0, 1); % the function value of the archived solutions

                  %% main loop
                  while nfes < max_nfes
                        pop = popold; % the old population becomes the current population
                        [fitness, sorted_index] = sort(fitness, 'ascend');
                        pop=pop(sorted_index,:);

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

                        r0 = 1 : pop_size;                
                        popAll = [pop; archive.pop];
                        [r1, r2, r3] = gnR1R2(pop_size, size(popAll, 1), r0);
                        vi = zeros(pop_size, problem_size);

                        %% mutation
                        bb = rand(pop_size, 1);
                        probiter = probDE(1,:);
                        l2 = sum(probDE(1:2));

                        de_1 = bb <= probiter(1)*ones(pop_size, 1);
                        de_2 = bb > probiter(1)*ones(pop_size, 1) &  bb <= (l2*ones(pop_size, 1));
                        de_3 = bb > l2*ones(pop_size, 1) &  bb <= (ones(pop_size, 1));

                        pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
                        randindex = floor(rand(1, pop_size) .* pNP) + 1; %% select from [1, 2, 3, ..., pNP]
                        pbest = pop(randindex, :); %% randomly choose one of the top 100p% solutions

                        % DE/current-to-p-best/1
                        vi(de_1==1, :)  = pop(de_1==1, :) + ...
                              sf(de_1==1, ones(1, problem_size)) .* ...
                              (pbest(de_1==1, :) - pop(de_1==1, :) + pop(r1(de_1==1), :) - popAll(r2(de_1==1), :));

                        % DE/rand/1 
                        vi(de_2==1, :) = pop(de_2==1, :) + ...
                              sf(de_2==1, ones(1, problem_size)) .* ...
                              (pop(r1(de_2==1), :) - popAll(r2(de_2==1), :));

                        % DE/current-to-q-best/1 with attraction
                        q_best_rate = 2 * p_best_rate - p_best_rate * (nfes/max_nfes);
                        qNP = max(round(q_best_rate * pop_size), 2); %% choose at least two best solutions
                        randindex = floor(rand(1, pop_size) .* qNP) + 1; %% select from [1, 2, 3, ..., pNP]
                        qbest = pop(randindex, :); %% randomly choose one of the top 100p% solutions

                        attraction = repmat(0.5 + 0.5 * (nfes/max_nfes), pop_size, problem_size); 

                        vi(de_3==1, :) = sf(de_3==1, ones(1, problem_size)) .* ...
                              (pop(r1(de_3==1), :) + ...
                              attraction(de_3==1, :) .* (qbest(de_3==1, :) - pop(r3(de_3==1), :)));

                        vi = boundConstraint(vi, pop, lu);

                        %% Doing q-best binomial crossover
                        mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
                        rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
                        jrand = sub2ind([pop_size problem_size], rows, cols);
                        mask(jrand) = false;

                        qNP = max(round(q_best_rate * size(popAll, 1)), 2); %% choose at least two best solutions
                        randindex = floor(rand(1, size(popAll, 1)) .* qNP) + 1; %% select from [1, 2, 3, ..., pNP]
                        popAllbest = popAll(randindex, :); %% randomly choose one of the top 100p% solutions
                        popAllbest = popAllbest(1:pop_size, :);

                        bb = rand(pop_size, 1) <= repmat(q_cr_rate, pop_size, 1);
                        qbest = pop;  qbest(bb, :) = popAllbest(bb, :);

                        ui = vi;      ui(mask) = qbest(mask);

                        children_fitness = feval(fhd, ui', func_no, C(cc,:));
                        children_fitness = children_fitness';

                          for i = 1 : pop_size
                              nfes = nfes + 1;

                              if children_fitness(i) < bsf_fit_var
                                  bsf_fit_var = children_fitness(i);
                                  bsf_solution = ui(i, :);                    
                              end

                              if nfes == max_nfes
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
                          %% ==================== update Prob. of each DE ===========================
                          diff2 = max(0,(fitness - children_fitness))./abs(fitness);
                          count_S(1)=max(0,mean(diff2(de_1==1)));
                          count_S(2)=max(0,mean(diff2(de_2==1)));
                          count_S(3)=max(0,mean(diff2(de_3==1)));

                          %% update probs.
                          if count_S~=0
                                probDE = max(0.1,min(0.9,count_S./(sum(count_S))));
                          else
                                probDE = 1.0/3 * ones(1,3);
                          end
                          %% ==================== update population and fitness =======================
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
                          else
                              memory_sf(memory_pos) = 0.5;
                              memory_cr(memory_pos) = 0.5;

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
                                  archive.funvalues = archive.funvalues(rndpos, :);           
                                end
                          end
                  end

                  if(C(cc,1)==1)
                        run_funcvals=run_funcvals-optimum(func_no);
                  end

                  run_funcvals(run_funcvals<val_2_reach)=0;
                  allerrorvals(:, run_id) = run_funcvals;
                  allsolutions(:, run_id) = bsf_solution';

            end %% end 1 run
            
            [~, sorted_index] = sort(allerrorvals(end,:), 'ascend');
            allerrorvals = allerrorvals(:, sorted_index);
            allsolutions = allsolutions(:, sorted_index);
            % Saving the solutions
            file_name=sprintf('Results/%s_%s_%s.txt',Alg_Name,int2str(func_no),int2str(problem_size));
            save(file_name, 'allerrorvals', '-ascii');
            file_name=sprintf('Solutions/%s_%s_%s.txt',Alg_Name,int2str(func_no),int2str(problem_size));
            save(file_name, 'allsolutions', '-ascii');
        end %% end 1 function run
    end %% end of a transformation & all functions
end %% end of the dimension
delete(gcp('nocreate'));   % close parallel pool to prevent memory overflow

end
