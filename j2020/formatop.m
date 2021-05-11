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


%% Read all datapoints & calculate required score data
for D=[10, 20]
    for f=1:10
        for i=1:size(C,1)
            Alg_Name  =['j2020_(' num2str(C(i,1)) num2str(C(i,2)) num2str(C(i,3)) ')'];
            file_path_C =['Results_C/' Alg_Name '_' num2str(f) '_' num2str(D) '.txt'];
            data=load(file_path_C);
            file_path_M =['Results/' Alg_Name '_' num2str(f) '_' num2str(D) '.txt'];
            save(file_path_M, 'data', '-ascii');
        end
    end
end