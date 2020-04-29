function p_value = one_sample_t_test(Data,condition,mu)

    if nargin < 3
        mu = 0;
        if nargin < 2
            condition = '==';
            if nargin < 1
                error 'Must give inputs...'
            end
        end
    end
    
    % Get T value
    mean_data = mean(Data);
    std_data = std(Data);
    n = length(Data);
    
    T_value = (mean_data-mu)/(std_data/sqrt(n));
    
    switch condition
        case '=='
            if T_value < 0
                p_value = 2*tcdf(T_value,n-1);
            else
                p_value = 2*tcdf(-T_value,n-1);
            end
            
        case '<'
            p_value = tcdf(-T_value,n-1);
        case '>'
            p_value = tcdf(T_value,n-1);
    end
    
end
            
            