function fitness = Parametrized_benchmark_func(pop,fun_n,C)
%
% Parametrized_benchmark_func evaluates (pop) Population on function (fun_n) using CEC2021 benchmark
% According to  the Parametrized Selection Vector (C), where:
% C1 is Bias indicator (0=0 ,1=F*)
% C2 Shift indicator (0=0, 1=oi)
% C3 Rotation indicator (0=I 1=M)
% Call:
%    fitness = Parametrized_benchmark_func(pop,fun_n,C)
%
% Version: 0.1  Date: 2021/11/17
% Written by Anas A. Hadi (anas1401@gmail.com) and Abhishek Kumar (abhishek.kumar.eee13@iitbhu.ac.in)



switch num2str(C)
    
    case '0  0  0'
        fitness=cec21_basic_func(pop,fun_n);
    case '1  0  0'
        fitness=cec21_bias_func(pop,fun_n);
    case '0  1  0'
        fitness=cec21_shift_func(pop,fun_n);
    case '0  0  1'
        fitness=cec21_rot_func(pop,fun_n);
    case '1  1  0'
        fitness=cec21_bias_shift_func(pop,fun_n);
    case '1  0  1'
        fitness=cec21_bias_rot_func(pop,fun_n);
    case '0  1  1'
        fitness=cec21_shift_rot_func(pop,fun_n) ;
    case '1  1  1'
        fitness=cec21_bias_shift_rot_func(pop,fun_n);
    otherwise
        disp('Undefined Selection Vector')
        
end



end