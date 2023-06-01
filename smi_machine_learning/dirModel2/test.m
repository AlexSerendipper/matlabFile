p = [2 2 2 2 0 0 0 0 0 -1 -1 -1 -1 1 1 1 1 0 0 0 0 0 -1 -1 -1 -1];

N = length(p);
i = 2;
while i<N
    if(p(i)==0)
        pre = p(i-1);  % 记录前值
        % j = i + 1;
        for j = i+1:N
            if(p(j)~=0)
                break
            end
        end
        n = floor((j-i)/2);
        
        
        % 将0的前半部分设为前值
        for k = i:i+n-1
            p(k) = pre;
        end
        % 将0的后半部分设为前值取反
        for k = i+n : j-1
            p(k) = -pre;
        end
        i = i + n -1;
    end
    i = i+1;
end