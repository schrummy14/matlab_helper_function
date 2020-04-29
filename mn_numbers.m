function nums = mn_numbers(m,n,N)

    nums = zeros(N,1);
    nums(1) = m;
    nums(2) = n;
    
    for k = 3:N
        nums(k) = nums(k-2) + nums(k-1);
    end
    
end