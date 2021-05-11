clc;
clear all;

format long;
format compact;

Results_10=[];
Results_30=[];
Results_50=[];
Results_100=[];

% Method_Name = 'LSHADE_EpSin_V_1';
% Method_Name = 'LSHADE_without_Archive';
% Method_Name = 'Ali_New_ADE';
% Method_Name = 'Ali_New_ADE_Ver2';
% Method_Name = 'Ali_New_ADE_with_Archive';
% Method_Name = 'Ali_New_ADE_with_Archive_Ver02';
% Method_Name = 'Ali_New_ADE_Success_Children_Ver_01';
% Method_Name = 'Ali_New_ADE_Success_Children_Ver_02';
% Method_Name = 'Ali_New_ADE_Success_Children_Ver_03';
% Method_Name = 'Ali_New_ADE_Success_Children_Ver_05';
% Method_Name = 'Ali_New_ADE_with_Archive_Ver02_and_Success_Children_Ver_05';
% Method_Name = 'Ali_New_ADE_with_Archive_and_Success_Children_Ver_02';
% Method_Name = 'Ali_New_ADE_with_Archive_and_Success_Children_Ver_07';
% Method_Name = 'Ali_New_ADE_with_Archive_and_Success_Children_Ver_02_R2';
% Method_Name = 'Ali_New_ADE_Success_Children_Ver_05_R2';
% Method_Name = 'Ali_New_ADE_Ver2_H2_Pop10';
% Method_Name = 'LSHADE_and_ANDE_V_2_R1';
Method_Name = 'LSHADE_SPA';
% Method_Name = 'LSHADE_EpSin_V_1';
% Method_Name = 'AliADENewWithArcAndSucc_EpSin_NoGRW_V0';
% Method_Name = 'AliADENewWithArc_EpSin_NoGRW_V1';
% Method_Name = 'AliADENewWithArcAndSuccHist_EpSin_NoGRW_V1';
% Method_Name = 'AliADENewWithArcAndSucc_EpSin_V1';
% Method_Name = 'AliADENewWithArc_EpSin_V1';
% Method_Name = 'AAliADENewWithArcAndSuccHist_EpSin_NoGRW_V2_Old_Mut_V2';
% Method_Name = 'SH_LSHADE_V2';
% Method_Name = 'LSHADE_NewF_V0_R2';
% Method_Name = 'LSHADE_NewF_GRW_V0_R2';
% Method_Name = 'LSHADE_NewF_LS_V0';
% Method_Name = 'LSHADE_NewF_LS_V1';














for func = 28:30
%     problem_size =10;
%     file_name=sprintf('Results\\%s_CEC2017_Problem#%s_problem_size#%s',Method_Name,int2str(func),int2str(problem_size));
%     load(file_name);
%     Results_10(func,:)=[min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome)];
% 
%     problem_size =30;
%     file_name=sprintf('Results\\%s_CEC2017_Problem#%s_problem_size#%s',Method_Name,int2str(func),int2str(problem_size));
%     load(file_name);
%     Results_30(func,:)=[min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome)];
%     
%     problem_size =50;
%     file_name=sprintf('Results\\%s_CEC2017_Problem#%s_problem_size#%s',Method_Name,int2str(func),int2str(problem_size));
%     load(file_name);
%     Results_50(func,:)=[min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome)];
    
    problem_size =100;
    file_name=sprintf('Results\\%s_CEC2017_Problem#%s_problem_size#%s',Method_Name,int2str(func),int2str(problem_size));
    load(file_name);
    Results_100(func,:)=[min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome)];
    
end %% end 1 function run

All_Results=[Results_10,Results_30,Results_50,Results_100];
All_Results(All_Results<10^-8)=0;
    

