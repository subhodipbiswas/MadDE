function [R1, R2, R3] = Gained_Shared_Junior_R1R2R3(indBest)


pop_size = length(indBest);
R0=1:pop_size;
R1=[];
R2=[];
R3=[];

for i=1:pop_size
    ind=find(indBest==i);
    if(ind==1)% best
    R1(i)=indBest(2);
    R2(i)=indBest(3);
    elseif(ind==pop_size)% worst
    R1(i)=indBest(pop_size-2);
    R2(i)=indBest(pop_size-1);
    else
    R1(i)=indBest(ind-1);
    R2(i)=indBest(ind+1);
    end
end

R3 = floor(rand(1, pop_size) * pop_size) + 1;

for i = 1 : 99999999
    pos = ((R3 == R2) | (R3 == R1) | (R3 == R0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
        R3(pos) = floor(rand(1, sum(pos)) * pop_size) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate R3 in 1000 iterations');
    end
end

end
