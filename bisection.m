function [h,fval,its] = bisection(fun,a,b,tol,maxIts,DOPLOT)

    if nargin < 6
        DOPLOT = false;
        if nargin < 5
            maxIts = 500;
            if nargin < 4
                tol = eps;
            end
        end
    end

    if nargin == 0
        fun = @(x) (x < 1e-7) - (x >= 1e-7);
        a = 1e-10;
        b = 1e-5;
        tol = eps;
        maxIts = 500;
        DOPLOT = true;
    end
    
    fa = fun(a);
    fb = fun(b);
    if fa*fb > 0
        error('Need a and b to be on opposite sides of the zero')
    end
    
    if DOPLOT
        textNum = @(k)['x',num2str(k)];
        offsetx = 0.1;
        offsety = 0.1;
        semilogx(logspace(-10,-5,1001),fun(logspace(-10,-5,1001)));
        xlim([1e-11, 1e-4])
        ylim(1.5*[-1 1])
        hold on
        plot(a,fun(a),'ko',b,fun(b),'ko')
        text(a*(1-offsetx),fun(a)*(1+offsety),'xa');
        text(b*(1-offsetx),fun(b)*(1-offsety),'xb');
        maxXs = 5;
    end
    for its = 1:maxIts
%         h = 0.5*(a+b);
        h = logspace(log10(a),log10(b),3);
        h = h(2);
        fval = fun(h);
        if fval*fb < 0
            a = h;
            fa = fval;
            if DOPLOT && its < maxXs
                plot(a,fa,'ko')
                text(a*(1-offsetx),fa*(1+offsety),textNum(its));
            end
        elseif fa*fval < 0
            b = h;
            fb = fval;
            if DOPLOT && its < maxXs
                plot(b,fb,'ko')
                text(b*(1-offsetx),fb*(1-offsety),textNum(its));
            end
        else
            break
        end
        if abs(a-b) < tol
            break
        end
    end
    if DOPLOT
        xlabel 'DEM Time Step (Seconds)'
        ylabel 'Function Output'
        hold off
    end
end

    
    