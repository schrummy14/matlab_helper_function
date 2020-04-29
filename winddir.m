function [w,j] = winddir(u,v,j) % ,method)

w = zeros(length(u),1);

switch nargin
    case 1 
        if isreal(u)
            error('Input must have two vectors or a single complex vector')
        end
        v = imag(u);
        u = real(u);
        w = winddir(u,v);
    case 2
        
        if length(v) == 1
            if isreal(u)
                error('Input must have two vectors or a single complex vector')
            end
            j = v;
            v = imag(u);
            u = real(u);
            [w,j] = winddir(u,v,j);
        else
            if numel(u)-numel(v)~=0
                error('vectors must be the same length')
            end
            for i = 1:length(u)
                if( u(i) >= 0 && v(i) > 0)
                    w(i) = mod(180/pi*abs(atan(u(i)/v(i))) + 270 + 180,360);
                elseif( u(i) < 0 && v(i) >= 0)
                    w(i) = mod(180/pi*abs(atan(v(i)/u(i))) + 180 + 180,360);
                elseif( u(i) <= 0 && v(i) < 0)
                    w(i) = mod(180/pi*abs(atan(u(i)/v(i))) +  90 + 180,360);
                elseif( u(i) > 0 && v(i) <= 0)
                    w(i) = mod(180/pi*abs(atan(v(i)/u(i))) +   0 + 180,360);
                end
            end
        end
        
    case 3
        
        if numel(u)-numel(v)~=0
            error('vectors must be the same length')
        end

        for i = 1:length(u)
    
            if( u(i) >= 0 && v(i) > 0)
                w(i) = mod(180/pi*abs(atan(u(i)/v(i))) + 270 + 180,360);
            elseif( u(i) < 0 && v(i) >= 0)
                w(i) = mod(180/pi*abs(atan(v(i)/u(i))) + 180 + 180,360);
            elseif( u(i) <= 0 && v(i) < 0)
                w(i) = mod(180/pi*abs(atan(u(i)/v(i))) +  90 + 180,360);
            elseif( u(i) > 0 && v(i) <= 0)
                w(i) = mod(180/pi*abs(atan(v(i)/u(i))) +   0 + 180,360);
            end
    
            w(i) = w(i) + j*360;
 
                if i > 1
                    if mod(w(i-1),360) > 330 && mod(w(i),360) < 30
                        w(i) = w(i) + 360;
                        j = j + 1;
                    elseif mod(w(i-1),360) < 30 && mod(w(i),360) > 330
                        w(i) = w(i) - 360;
                        j = j - 1;
                    end 
                end
        end
        
%     case 4
%         
%         if numel(u)-numel(v)~=0
%             error('vectors must be the same length')
%         end
% 
%         for i = 1:length(u)
%     
%             if( u(i) >= 0 && v(i) > 0)
%                 w(i) = mod(180/pi*abs(atan(u(i)/v(i))) + 270 + 180,360);
%             elseif( u(i) < 0 && v(i) >= 0)
%                 w(i) = mod(180/pi*abs(atan(v(i)/u(i))) + 180 + 180,360);
%             elseif( u(i) <= 0 && v(i) < 0)
%                 w(i) = mod(180/pi*abs(atan(u(i)/v(i))) +  90 + 180,360);
%             elseif( u(i) > 0 && v(i) <= 0)
%                 w(i) = mod(180/pi*abs(atan(v(i)/u(i))) +   0 + 180,360);
%             end
%     
%             w(i) = w(i) + j*360;
%             
%             if method == 1
%                 if i > 1
%                     a = w(i)+360 - w(i-1);
%                     b = w(i)-360 - w(i-1);
%                     c = w(i)-w(i-1);
% 
%                     if abs(a)<abs(b) && abs(a)<abs(c)
%                         w(i) = w(i)+360;
%                         j = j + 1;
%                     elseif abs(b)<abs(a) && abs(b)<abs(c)
%                         w(i) = w(i)-360;
%                         j = j - 1;
%                     end
%                 end
%             else
%                 if i > 1
%                     if mod(w(i-1),360) > 330 && mod(w(i),360) < 30
%                         w(i) = w(i) + 360;
%                         j = j + 1;
%                     elseif mod(w(i-1),360) < 30 && mod(w(i),360) > 330
%                         w(i) = w(i) - 360;
%                         j = j - 1;
%                     end 
%                 end
%             end
%         end
end
        