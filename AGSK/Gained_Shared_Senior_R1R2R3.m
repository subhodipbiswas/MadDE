function [R1, R2, R3] = Gained_Shared_Senior_R1R2R3(indBest)


pop_size = length(indBest);

R1=indBest(1:round(pop_size*0.05));
R1rand = ceil(length(R1) * rand(pop_size, 1));
R1 = R1(R1rand);

R2=indBest(round(pop_size*0.05)+1:round(pop_size*0.95));
R2rand = ceil(length(R2) * rand(pop_size, 1));
R2 = R2(R2rand);

R3=indBest(round(pop_size*0.95)+1:end);
R3rand = ceil(length(R3) * rand(pop_size, 1));
R3 = R3(R3rand);

end
