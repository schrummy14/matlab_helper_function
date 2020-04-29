function a = mag(x,y,z,varargin)

switch nargin
    case 1
        if size(x,2) > 1
            newVar = num2cell(x,1);
            a = mag(newVar{:});
        else
            a=abs(x);
        end
    case 2
        if numel(x) == numel(y)
            a=sqrt(x.^2+y.^2);
        else
            error('Check sizes of entries...')
        end
    case 3
        if numel(x) == numel(y) && numel(x) == numel(z)
            a=sqrt(x.^2+y.^2+z.^2);
        else
            error('Check sizes of entries...')
        end
    otherwise
        b = zeros(numel(varargin),1);
        for n = 1:numel(varargin)
            b(n) = numel(varargin{n});
        end
        c = sum(b)/length(b);
        
        if numel(x) == numel(y) && numel(x) == numel(z) && numel(x)==c 
            a1 = x.^2+y.^2+z.^2;
            a2 = 0;
            for n = 1:length(varargin)
                a2 = a2 + varargin{n}.^2;
            end
            a = sqrt(a1 + a2);
        else
            error('Check sizes of entries...')
        end
end
