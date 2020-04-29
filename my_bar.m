function [h_b, h_l] = my_bar(data,numBars,xlabel_str,ylabel_str,title_str)

    if nargin < 5
        title_str = [];
        if nargin < 4
            ylabel_str = 'Percentage of Points';
            if nargin < 3
                xlabel_str = 'Observed Values';
                if nargin < 2
                    numBars = 9;
                    if nargin < 1
                        error 'You must at least supply the data...'
                    end
                end
            end
        end
    end
    
%     [a,b] = hist(data,numBars);
%     
%     [m,n] = size(a);
%     if min(m,n) == 1
%         multi_obs = false;
%         h_b = bar(b,100*a/sum(a));
%     else
%         multi_obs = true;
%         h_b = bar(b,100*a/sum(a(:,1)));
%     end
%     for k = 1:length(h_b)
%         h_b(k).Parent.XTick = b;
%     end
%     
%     title(title_str);
%     xlabel(xlabel_str);
%     ylabel(ylabel_str);
%     
%     if multi_obs
%         h_l = legend;
%     else
%         h_l = nan;
%     end
    
    h_b = histogram(data,numBars);
    h_b.Normalization = 'probability';
    hValues = h_b.Parent.YTick;
    hYTickNames = cell(length(hValues),1);
    for k = 1:length(hValues)
        hYTickNames{k} = num2str(100*hValues(k));
    end
    h_b.Parent.YTickLabel = hYTickNames;
    h_b.BinWidth = 0.95*h_b.BinWidth;
    title(title_str);
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    
    
    h_l = nan;

end