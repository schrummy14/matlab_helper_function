function dx = numDer(x,h,derNumber,order,type)
    if size(x,2)>1
        x = x.';
    end
    numX = length(x);
    switch type
        case 'central'
            switch derNumber
                case 1
                    if any(order == [1,3,5,7])
                        error('central methods are at least order 2n...')
                    end
                    if order > 8
                        warning('Order is greater than max order method, using 8th order...')
                        order = 8;
                    end
                    dxCoef = [
                        0       0       0       0       0       0       0       0       0
                        0       0       0      -1/2     0       1/2     0       0       0
                        0       0       0       0       0       0       0       0       0
                        0       0       1/12   -2/3     0       2/3    -1/12    0       0
                        0       0       0       0       0       0       0       0       0
                        0      -1/60    3/20   -3/4     0       3/4    -3/20    1/60    0
                        0       0       0       0       0       0       0       0       0
                        1/280  -4/105   1/5    -4/5     0       4/5    -1/5     4/105  -1/280];
                    
                case 2
                    if any(order == [1,3,5,7])
                        error('central methods are at least order 2n...')
                    end
                    if order > 8
                        warning('Order is greater than max order method, using 8th order...')
                        order = 8;
                    end
                    
                    dxCoef = [
                        0       0       0       0       0       0       0       0       0
                        0       0       0       1      -2       1       0       0       0
                        0       0       0       0       0       0       0       0       0
                        0       0      -1/12    4/3    -5/2     4/3    -1/12    0       0
                        0       0       0       0       0       0       0       0       0
                        0       1/90   -3/20    3/2    -49/18   3/2    -3/20    1/90    0
                        0       0       0       0       0       0       0       0       0
                       -1/560   8/315  -1/5     8/5    -205/72  8/5    -1/5     8/315  -1/560];
                   
            end
            
            if length(x)~=9
                g = zeros(9,1);
                glow = 5-ceil(length(x)/2)+1;
                ghig = 5+ceil(length(x)/2)-1;
                g(glow:ghig) = x;
                x = g;
                if numX<order+1
                    warning('Not enough points to use order... using max order...')
                    order = numX-1;
                end
            end
            
        case 'forward'
            switch derNumber
                case 1
                    dxCoef = [
                       -1       1       0       0       0       0       0   0
                       -3/2     2      -1/2     0       0       0       0   0
                       -11/6    3      -3/2     1/3     0       0       0   0
                       -25/12   4      -3       4/3    -1/4     0       0   0
                       -137/60  5      -5       10/3   -5/4     1/5     0   0
                       -49/20   6      -15/2    20/3   -15/4    6/5    -1/6 0
                       ];
                   if order > 6
                       warning('Order is greater than max order method, using 6th order...')
                       order = 6;
                   end
                case 2
                    dxCoef = [
                        1         -2       1         0          0        0          0           0
                        2         -5       4        -1          0        0          0           0
                        35/12     -26/3    19/2     -14/3       11/12    0          0           0
                        15/4      -77/6    107/6    -13         61/12   -5/6        0           0
                        203/45    -87/5    117/4    -254/9      33/2    -27/5       137/180     0
                        469/90    -223/10  879/20   -949/18     41      -201/10     1019/180   -7/10];
                    
            end
           
           
        case 'backward'
            dxCoef = [
                0   0       0       0       0       0      -1       1
                0   0       0       0       0       1/2    -2       3/2
                0   0       0       0      -1/3     3/2    -3       11/6
                0   0       0       1/4    -4/3     3      -4       25/12
                0   0      -1/5     5/4    -10/3    5      -5       137/60
                0   1/6    -6/5     15/4   -20/3    15/2   -6       49/20   ];
            if order > 6
               warning('Order is greater than max order method, using 6th order...')
               order = 6;
            end
           
    end
    
    dx = dxCoef(order,:)*x/h^derNumber;

end


