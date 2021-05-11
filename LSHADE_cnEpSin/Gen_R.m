function R = Gen_R(NP_Size,N)

% Gen_R generate N column vectors r1, r2, ..., rN of size NP_Size
%    R's elements are choosen from {1, 2, ..., NP_Size} & R(j,i) are unique per row

% Call:
%    [R] = Gen_R(NP_Size)   % N is set to be 1;
%    [R] = Gen_R(NP_Size,N) 
%
% Version: 0.1  Date: 2018/02/01
% Written by Anas A. Hadi (anas1401@gmail.com)


R(1,:)=1:NP_Size;

for i=2:N+1
    
    R(i,:) = ceil(rand(NP_Size,1) * NP_Size);
    
    flag=0;
    while flag ~= 1
        pos = (R(i,:) == R(1,:));
        for w=2:i-1
            pos=or(pos,(R(i,:) == R(w,:)));
        end
        if sum(pos) == 0
            flag=1;
        else
            R(i,pos)= floor(rand(sum(pos),1 ) * NP_Size) + 1;
        end
    end
end

R=R';

end
