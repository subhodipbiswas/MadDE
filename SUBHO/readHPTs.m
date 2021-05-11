function [parm_settings] = readHPTs(folder_name)
% Read hyperparameters from existing configuration files. These files are
% stored in the folder 'Hyperparameters'.

ctr = 0;
parm_settings = zeros(ctr,9);

getFolderInfo = dir(folder_name);
[num_files, ~] = size(getFolderInfo);

if num_files > 2
    % There are actual stored values of hyperparameters present from a
    % previous run of SUBHO
    for file_num = 3:num_files
        % f_path = [getFolderInfo(file_num).folder '/' getFolderInfo(file_num).name ];
        f_path = [getFolderInfo(file_num).folder '\' getFolderInfo(file_num).name ];
        if isfile(f_path)
            fprintf('\nReading: %s', f_path);
            ctr = ctr+1;
            parm_settings(ctr,:) = load(f_path);
        end
    end
    fprintf('\n\nUsing existing hyperparameters for warm-start...\n\n')
else
    fprintf('\n\nDid not find existing hyperparameters. Start from scratch...\n\n')
end

end
