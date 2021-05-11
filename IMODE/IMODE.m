%% ============ Improved Multi-operator Differential Evolution Algorithm (IMODE) ============
% Should you have any queries, please contact
% Dr. Karam Sallam. Zagazig University
% karam_sallam@zu.edu.eg
% =========================================================================
function [x, xold, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq,current_eval,res_det,F,cr ] = ...
    IMODE( x,xold, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq, xmin, xmax,n,...
    PopSize,current_eval,I_fno,C,res_det,Printing,Max_FES, G_Max, gg,F,cr)

vi=zeros(PopSize,n);

%% calc CR and F
mem_rand_index = ceil(memory_size * rand(PopSize, 1));
mu_sf = archive_f(mem_rand_index);
mu_cr = archive_Cr(mem_rand_index);

%% ========================= generate CR ==================================

cr = normrnd(mu_cr, 0.1);
term_pos = find(mu_cr == -1);
cr(term_pos) = 0;
cr = min(cr, 1);
cr = max(cr, 0);
%         cr=0.1.*ones(1,PopSize);
%% for generating scaling factor


F = mu_sf + 0.1 * tan(pi * (rand(1,PopSize) - 0.5));
pos = find(F <= 0);

while ~ isempty(pos)
    F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(1,length(pos)) - 0.5));
    pos = find(F <= 0);
end

F = min(F, 1);
F=F';
[fitx,inddd]=sort(fitx);
x=x(inddd,:);
[cr,~]=sort(cr);


%% ======================== generate new x =================================
popAll = [x;archive.pop]; %% set archive
r0 = 1 : PopSize;
%% generate random integer numbers
[r1, r2,r3] = gnR1R2(PopSize, size(popAll, 1), r0);

%% mutation
bb= rand(PopSize, 1);
probiter = prob(1,:);
l2= sum(prob(1:2));
op_1 = bb <=  probiter(1)*ones(PopSize, 1);
op_2 = bb > probiter(1)*ones(PopSize, 1) &  bb <= (l2*ones(PopSize, 1)) ;
op_3 = bb > l2*ones(PopSize, 1) &  bb <= (ones(PopSize, 1)) ;

[~,inddd]=sort(fitx);

pNP = max(round(0.25* PopSize), 1); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
phix = x(randindex, :);
vi(op_1==1,:) = x(op_1==1,:)+ F(op_1==1, ones(1, n)) .*(phix(op_1==1,:) - x(op_1==1,:) + x(r1(op_1==1), :) - popAll(r2(op_1==1), :));
vi(op_2==1,:) =  x(op_2==1,:)+ F(op_2==1, ones(1, n)) .*(phix(op_2==1,:) - x(op_2==1,:) + x(r1(op_2==1), :) - x(r3(op_2==1), :));
%% DE3
pNP = max(round(0.5 * PopSize), 2); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
phix = x(randindex, :);
vi(op_3==1,:) = F(op_3==1, ones(1, n)).* x(r1(op_3==1), :) +F(op_3==1, ones(1, n)).* (phix(op_3==1,:) - x(r3(op_3==1), :));

%% handle boundaries
vi = han_boun(vi, xmax, xmin, x,PopSize,2);
%% crossover
if rand<0.4
    mask = rand(PopSize, n) > cr(:, ones(1, n));
    % mask is used to indicate which elements of ui comes from the parent
    rows = (1 : PopSize)'; cols = floor(rand(PopSize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([PopSize n], rows, cols); mask(jrand) = false;
    ui = vi; ui(mask) = x(mask);
else
    ui=x;
    startLoc= randi(n,PopSize,1);
    for i=1:PopSize
        l=startLoc(i);
        while (rand<cr(i) && l<n)
            l=l+1;
        end
        for j=startLoc(i) : l
            ui(i,j)= vi(i,j);
        end
    end
end
% ui = x; ui(mask) = vi(mask);
%% evaluate
fitx_new = Parametrized_benchmark_func(ui',I_fno,C);
%% update FITNESS EVALUATIONS
current_eval =current_eval+PopSize;

%% calc. imprv. for Cr and F
diff = abs(fitx - fitx_new);
I =(fitx_new < fitx);
goodCR = cr(I == 1);
goodF = F(I == 1);

%% ========================= update archive ===============================
archive = updateArchive(archive, x(I == 1, :), fitx(I == 1)');
%% ==================== update Prob. of each DE ===========================
diff2 = max(0,(fitx - fitx_new))./abs(fitx);
count_S(1)=max(0,mean(diff2(op_1==1)));
count_S(2)=max(0,mean(diff2(op_2==1)));
count_S(3)=max(0,mean(diff2(op_3==1)));

%% update probs.
if count_S~=0
    prob= max(0.1,min(0.9,count_S./(sum(count_S))));
else
    prob=1/3 * ones(1,3);
end
%% ==================== update x and fitx =================================
fitx(I==1)= fitx_new(I==1);
xold(I == 1, :) = x(I == 1, :);
x(I == 1, :) = ui(I == 1, :);

%% =================== update memory cr and F =============================

if size(goodF,1)==1
    goodF=goodF';
end
if size(goodCR,1)==1
    goodCR=goodCR';
end
num_success_params = numel(goodCR);
if num_success_params > 0
    weightsDE = diff(I == 1)./ sum(diff(I == 1));
    %% for updating the memory of scaling factor
    archive_f(hist_pos) = (weightsDE * (goodF .^ 2))./ (weightsDE * goodF);
    
    %% for updating the memory of crossover rate
    if max(goodCR) == 0 || archive_Cr(hist_pos)  == -1
        archive_Cr(hist_pos)  = -1;
    else
        archive_Cr(hist_pos) = (weightsDE * (goodCR .^ 2)) / (weightsDE * goodCR);
    end
    
    hist_pos= hist_pos+1;
    if hist_pos > memory_size;  hist_pos = 1; end
else
    archive_Cr(hist_pos)=0.5;
    archive_f(hist_pos)=0.5;
    % end
end

%% sort new x, fitness
[fitx, ind]=sort(fitx);
x=x(ind,:);
xold = xold(ind,:);

%% record the best value after checking its feasiblity status
if fitx(1)<bestold  && min(x(ind(1),:))>=-100 && max(x(ind(1),:))<=100
    bestold=fitx(1);
    bestx= x(1,:);
end
%% check to print
if Printing==1
    res_det= [res_det repmat(bestold,1,PopSize)];
end

