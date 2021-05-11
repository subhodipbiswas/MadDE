%% Script to calculate Scores on CEC 2021 Benchmark Suite
%  ************************************************************************
%  For any queries please feel free to contact us
%  Subhodip Biswas : sub17was(at)gmail.com
%  Debanjan Saha   : debanjansh(at)gmail.com
%  %%%%%%%%%%% Version Log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Version         : 1.0.0        Packaged as part of Initial Code Release                
%  Date-Written    : 2021/01/23   of MadDE at CEC2021 SO-BCO       
%  ************************************************************************

clc;
clear;
fprintf('\n\n******* SCORES  ON  CEC2021 single-objective optimization ******\n\n Algorithms:\t| ');
% Add/Modify the algo to include all the comparing algos
algo={
    ['AGSK'],...
    ['IMODE'],...
    ['j2020'],...
    ['LSHADE'],...
    ['LSHADE_cnEpSin'],...
    ...['LSHADE_RSP'],...
    ['MadDE']
    };
n_algo=size(algo,2); % no of algorithms being compared
for a=1:n_algo
    fprintf(" %s | ",char(algo(a)));
end
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
                file_path =[char(algo(j)) '/Results/' Alg_Name '_' num2str(f) '_' num2str(D(k)) '.txt'];
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
    score1(ind)=50*(min_sne/sne(ind));
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
    score2(ind)=50*(min_sr/sr(ind));
    ind = ind+1;
end

%% Display Final Score
print_res=zeros(n_algo,1);
for i=1:n_algo
    score(i)=(score1(i)+score2(i));
    fprintf("\n\nSNE for %s : %.2f\n",char(algo(i)),sne(i));
    fprintf("SR for %s : %.2f\n",char(algo(i)),sr(i));
    fprintf("Score 1 for %s : %.2f\n",char(algo(i)),score1(i));
    fprintf("Score 2 for %s : %.2f\n",char(algo(i)),score2(i));
    fprintf("\nOverall Score for %s : %.2f\n",char(algo(i)),score(i));
    
end
fprintf("*************************************************\n");
