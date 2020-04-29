function fizzbuzz(n)   

    if nargin<1
        n = 25;
    end

    num_check = [
        3
        5
        ];
    
    output_var= {
        'fizz'
        'buzz'
        };
    
    parfor m = 1:n
        output = ['(',num2str(m),') '];
        for k = 1:length(num_check)
            if mod(m,num_check(k))==0
                output = [output,output_var{k}]; 
            end
        end
        if strcmp(output,'')
            output = num2str(m);
        end
        disp(output)
    end

end