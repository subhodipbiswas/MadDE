function [score] = tune_hyperparameters(var, num_hyper)
%EVAL_OBJECTIVE Summary of this function goes here
%   Detailed explanation goes here


for i=1:num_hyper
    switch i
        case 1
            varname = 'q_cr_rate';
        case 2
            varname = 'p_best_rate';
        case 3
            varname = 'arc_rate';
        case 4
            varname = 'mem_mult';
        case 5
            varname = 'pop_mult';
        case 6
            varname = 'sf_init';
        case 7
            varname = 'cr_init';
        otherwise
            disp('Error..')
    end
    
    str = sprintf('%s= var.%s;', varname, varname);
    eval(str)
end

MadDE(q_cr_rate, p_best_rate, arc_rate, mem_mult, pop_mult,...
    sf_init, cr_init, 'MadDEvHPT');

% compute score of this parameteric configuration w.r.t. MadDEv0.0
[score_1, score_2] = getScoreopt();
score = score_2 + (100 - score_1);

% save files
save_res = zeros(9,1);
save_res(1,:) = q_cr_rate;
save_res(2,:) = p_best_rate;
save_res(3,:) = arc_rate;
save_res(4,:) = mem_mult;
save_res(5,:) = pop_mult;
save_res(6,:) = sf_init;
save_res(7,:) = cr_init;
save_res(8,:) = score_1;
save_res(9,:) = score_2;
out_format = 'ddmmyy_HH_MM';
out_fname = datestr(now,out_format);
file_name=sprintf('Hyperparameters/HPTSettings_%s.txt',out_fname);
save(file_name, 'save_res', '-ascii');

end
