function need_primes = get_primes(minNum, numPrimes)
    need_primes = zeros(numPrimes,1);
    k = 0;
    n = ceil((minNum-1)/6);
    while true
        maybe_rand = 6*n-1;
        if isprime(maybe_rand)
            k = k + 1;            
            need_primes(k,1) = maybe_rand;
            if k >= numPrimes
                break
            end
        end
        maybe_rand = 6*n+1;
        if isprime(maybe_rand)
            k = k + 1;            
            need_primes(k,1) = maybe_rand;
            if k >= numPrimes
                break
            end
        end
        n = n + 1;
    end
end