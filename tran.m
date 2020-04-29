function A = tran(A)
    n = length(A);
    for i=1:n
        for j=i:n
            if i~=j
                temp = A(i,j);
                A(i,j) = A(j,i);
                A(j,i) = temp;
%                 disp(A)
            end
        end
    end
end

            