function [r1, r2,r3,r4,r5,r6,r7] =generate_random(NP1)
r1=zeros(1,NP1);
r2=zeros(1,NP1);
r3=zeros(1,NP1);
r4=zeros(1,NP1);
r5=zeros(1,NP1);
r6=zeros(1,NP1);
r7=zeros(1,NP1);
for i=1:NP1
    r1(i)=ceil(rand*NP1);
    while  r1(i)==i
        r1(i)=ceil(rand*NP1);
    end
    r2(i)=ceil(rand*NP1);
    while  r2(i)==r1(i) || r2(i)==i
        r2(i)=ceil(rand*NP1);
    end
    r3(i)=ceil(rand*NP1);
    while r3(i)==r2(i) ||r3(i)==r1(i)|| r3(i)==i
        r3(i)=ceil(rand*NP1);
    end
    
    r4(i)=ceil(rand*NP1);
    while r4(i)==r2(i) ||r4(i)==r1(i)|| r4(i)==i ||r4(i)==r3(i)
        r4(i)=ceil(rand*NP1);
    end
    
    r5(i)=ceil(rand*NP1);
    while r5(i)==r2(i) ||r5(i)==r1(i)|| r5(i)==i ||r5(i)==r3(i) ||r5(i)==r4(i)
        r5(i)=ceil(rand*NP1);
    end
    
    
    
    r6(i)=ceil(rand*NP1);
    while r6(i)==r2(i) ||r6(i)==r1(i)|| r6(i)==i ||r6(i)==r3(i) ||r6(i)==r4 (i) ||r6(i)==r5(i)
        r6(i)=ceil(rand*NP1);
    end
    
      r7(i)=ceil(rand*NP1);
    while r7(i)==r2(i) ||r7(i)==r1(i)|| r7(i)==i ||r7(i)==r3(i) ||r7(i)==r4 (i) ||r7(i)==r5(i) ||r7(i)==r6(i)
        r7(i)=ceil(rand*NP1);
    end
    
    
end
end