%% ============ Improved Multi-operator Differential Evolution Algorithm (IMODE) ============
% Should you have any queries, please contact
% Dr. Karam Sallam. Zagazig University
% karam_sallam@zu.edu.eg
% =========================================================================
% Some part of this code is taken from UMOEA-II
% =========================================================================
function [Par] = Introd_Par(I_fno,D,C)

%% loading

Par.n_opr=3;  %% number of operators 
Par.n=D;     %% number of decision vriables


if Par.n==10
    Par.CS=100; %% cycle
    Par.Max_FES=200000;
    Par.Gmax = 2163;
elseif Par.n==20
    Par.CS=100; %% cycle
    Par.Gmax = 2745;
    Par.Max_FES=1000000;
else
    Par.Max_FES=10000000;
    Par.CS=100; %% cycle
    Par.Gmax = 3401;
end
optima= [100, 1100 ,700 ,1900 ,1700 ,1600 ,2100 ,2200 ,2400 ,2500]; %% define the optimal solution as shown in the TR
Par.xmin= -100*ones(1,Par.n);
Par.xmax= 100*ones(1,Par.n);

if C(1)==1
    optimum=optima(I_fno);
else
    optimum = 0;
end
Par.f_optimal=optimum;
Par.PopSize=6*Par.n*Par.n; %% population size
% Par.PopSize=18*Par.n; %% population size

Par.MinPopSize=4;
Par.MinPopSize1=4+floor(3*log(Par.n));

Par.prob_ls=0;    % For deactivating LS2
%% printing the detailed results- this will increase the computational time
Par.Printing=1; %% 1 to print; 0 otherwise

end