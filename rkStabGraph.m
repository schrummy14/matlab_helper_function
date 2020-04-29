function rkStabGraph(A,b,wind,mesh)

    if nargin == 0
%         A = [
%             0.0 0.0 0.0 0.0
%             0.5 0.0 0.0 0.0
%             0.0 0.5 0.0 0.0
%             0.0 0.0 1.0 0.0
%             ];
        A = [0 0 0
            5/24 1/3 -1/24
            1/6 2/3 1/6];
        b = [1/6 2/3 1/6];
    end

    x1 = -5;
    x2 =  5;
    y1 = -5;
    y2 =  5;
    nx = 101;
    ny = 101;

    if nargin > 2
        if numel(wind) > 2
            x1 = wind(1);
            x2 = wind(2);
            y1 = wind(3);
            y2 = wind(4);
        elseif numel(wind) > 1
            x2 = wind(1);
            x1 = -x2;
            y2 = wind(2);
            y1 = -y2;
        else
            x1 = -wind;
            x2 = wind;
            y1 = x1;
            y2 = x2;
        end
        if nargin > 3
            if numel(mesh) > 1
                nx = mesh(1);
                ny = mesh(2);
            else
                nx = mesh;
                ny = nx;
            end
        end
    end

    windx = linspace(x1,x2,nx);
    windy = linspace(y1,y2,ny);

    [X,Y] = meshgrid(windx,windy);

    Z = X+1i*Y;

    znew = reshape(Z,numel(Z),1);
    
    Rnew = zeros(nx,ny,size(b,1));
    for m = 1:size(b,1)
        
        s = size(b(m,:),2);

        e=ones(s,1);
        p=poly(A-e*b(m,:)); %Numerator
        p = p(end:-1:1);
        q=poly(A);         %Denominator
        q = q(end:-1:1);
        R1 = polyval(p,znew);
        R2 = polyval(q,znew);
        R = abs(R1./R2);

        Rnew(:,:,m) = reshape(R,nx,ny);
    end
    
    if size(b,1) < 2
        v = 0:.01:2;
        figure
        colormap(jet(length(v)))
        contourf(windx,windy,Rnew,v,'LineColor','None');
        caxis([0 2])
        hold on
        contour(windx,windy,Rnew,[1 1],'LineColor','Black','LineWidth',2);
        hold off
    else
        figure
        for m = 1:size(b,1)
            hold on
            contour(windx,windy,Rnew(:,:,m),[1 1],'LineColor','Black','LineWidth',2);
        end
        hold off
    end
end