function vi = boundConstraint (vi, pop,Par)
% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang, jingqiao@gmail.com

[NP, D] = size(pop);  % the population size and the problem's dimension
ht=randi(2); %% boundary handling type
if ht==1
%% check the lower bound
xl = repmat(Par.xmin.*ones(1,Par.n), NP, 1);
pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;

%% check the upper bound
xu = repmat(Par.xmax.*ones(1,Par.n), NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;
elseif ht==2
     x_L = repmat(Par.xmin.*ones(1,Par.n), NP, 1);
    x_U = repmat(Par.xmax.*ones(1,Par.n), NP, 1);
    pos = vi < x_L | vi>x_U;
    vi(pos) = x_L(pos) + rand*(x_U(pos) - x_L(pos));
    
end