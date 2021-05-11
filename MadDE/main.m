%% ************************************************************************
% This package is a MATLAB source code of CEC 2021 Bound Constrained Competition
% Execute proposed algorithm (MadDE) on CEC 2021 Benchmark
%% ************************************************************************
% About: CEC 2021 Bound Constrained Competition
% Ali Wagdy, Anas A Hadi, Ali K. Mohamed, Prachi Agrawal, Abhishek Kumar, and P. N. Suganthan,
% Problem Definitions and Evaluation Criteria for the CEC 2021 Special Session and Competition on Single Objective Bound Constrained Numerical Optimization,
% Technical Report (TR), Nanyang Technological University, Singapore.
% https://www3.ntu.edu.sg/home/epnsugan/index_files/CEC2021/CEC2021-2.htm
%
%% ************************************************************************
%  For any queries please feel free to contact us
%  Subhodip Biswas : sub17was(at)gmail.com
%  Debanjan Saha   : debanjansh(at)gmail.com
%  %%%%%%%%%%% Version Log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Version         : 1.0.0        Written code of MadDE                 
%  Date-Written    : 2021/01/23           
mex cec21_basic_func.cpp -DWINDOWS
mex cec21_bias_func.cpp -DWINDOWS
mex cec21_bias_rot_func.cpp -DWINDOWS
mex cec21_bias_shift_func.cpp -DWINDOWS
mex cec21_bias_shift_rot_func.cpp -DWINDOWS
mex cec21_rot_func.cpp -DWINDOWS
mex cec21_shift_func.cpp -DWINDOWS
mex cec21_shift_rot_func.cpp -DWINDOWS

clc;
clear all;

format long;
Runs=30;
fhd=@Parametrized_benchmark_func;

optimum= [100, 1100, 700, 1900, 1700, 1600, 2100, 2200, 2400, 2500]; % TR

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

func=[1:10]; %% 10 functions defined in the CEC 2021 TR
%%  Execute the algorithm per dimensions defined
for D= [10 20]
    switch D
        case 10
            max_nfes=200000;
        case 20
            max_nfes=1000000;
        otherwise
            disp('\nError!! Invalid Dimension Size Provided\n\n')
    end
    % for each transformation execute MadDE algorithm
    for i=1:size(C,1)
        fprintf('\n-------------------------------------------------------\n\n')
        MadDE(Runs,fhd,C(i,:),D,func,max_nfes,optimum);
    end
end

