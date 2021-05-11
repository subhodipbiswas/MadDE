function [score_1, score_2] = getScoreopt()
%COMPUTESCORE Summary of this function goes here
%   Returns difference in score for an optimal version & the Bayes version

pause(5);
algo={['MadDEvHPT'], ['MadDEv0.0']};
n_algo=size(algo,2); % no of algorithms being compared

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
% Initialize variables
xbest=zeros(160, n_algo);
xsr=zeros(160, n_algo);
D=[10, 20];

%% Read all datapoints & calculate required score data
for k=1:2
    for f=1:10
        for i=1:size(C,1)
            for j=1:n_algo
                % for each algo read the results
                Alg_Name  =[char(algo(j)) '_(' num2str(C(i,1)) num2str(C(i,2)) num2str(C(i,3)) ')'];
                file_path =['Results/' Alg_Name '_' num2str(f) '_' num2str(D(k)) '.txt'];
%                 fprintf("Trying to read file= %s\n",file_path);
                data=load(file_path);
                last_row=data(16,:);
                xbest((8*(f-1))+80*(k-1)+i,j)=min(last_row);
                xmean((8*(f-1))+80*(k-1)+i,j)=mean(last_row);
            end
            xbest_max((8*(f-1))+80*(k-1)+i,:)=max(xbest((8*(f-1))+80*(k-1)+i,:));
            xsr((8*(f-1))+80*(k-1)+i,:)=tiedrank(xmean((8*(f-1))+80*(k-1)+i,:));
            % perform normalization
            for j=1:n_algo
                if xbest((8*(f-1))+80*(k-1)+i,j) == 0
                    ne_arr((8*(f-1))+80*(k-1)+i,j) = 0;
                else
                    ne_arr((8*(f-1))+80*(k-1)+i,j) = xbest((8*(f-1))+80*(k-1)+i,j)/xbest_max((8*(f-1))+80*(k-1)+i,1);
                end
            end
        end
    end
end

%% Calculate SNE
%  SNE is defined as the average of all normalized error values over all 
%  functions, configurations and dimensions.
ind=1;
while ind <= n_algo
    sne(ind)=0.5*sum(ne_arr(:,ind));
    ind = ind+1;
end
min_sne = min(sne(1,:));

%% Calculate Score 1
%  It is the sum of normalized error values
ind=1;
while ind <= n_algo
    score1(ind)=100*(min_sne/sne(ind));
    ind = ind+1;
end

%% Calculate SR
%  SR is defined as the sum of ranks where rank is the algorithm�s rank 
%  among all algorithms for a given function, configuration and
%  dimension that is based on its mean error value (not normalized).
ind=1;
while ind <= n_algo
    sr(ind)=0.5*sum(xsr(:,ind));
    ind = ind+1;
end
min_sr = min(sr(1,:));

%% Calculate Score 2
%  It is the sum of ranks among all algorithms
ind=1;
while ind <= n_algo
    score2(ind)=100*(min_sr/sr(ind));
    ind = ind+1;
end

%% Display Final Score
print_res=zeros(n_algo,1);
for i=1:n_algo
    score(i) = 0.5 * (score1(i) + score2(i));
    print_res(i,1)=sne(i);
    print_res(i,2)=sr(i);
    print_res(i,3)=score1(i);
    print_res(i,4)=score2(i);
    print_res(i,5)=score(i);
end

% algo_score = score(2) - score(1);
score_1 = score(1);
score_2 = score(2);

end

