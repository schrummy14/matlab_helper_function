function [ varargout ] = subplot_position_manager( varargin )
%subplot_position_manager Gets number of subplots,form of them, and
%required margins to produce proper position of each subplot
LegPos = '';
LegSize = [0.2 0.2];
            HLeg = 0.2;
            VLeg = 0.2;
    for n = 1:2:numel(varargin)
        switch varargin{n}
            case 'PerRow'
                Row = varargin{n+1};
            case 'Margin'
                Marg = varargin{n+1};
            case 'WhiteSpace'
                Wh = varargin{n+1};
            case 'LegendPosition'
                LegPos = varargin{n+1};
            case 'LegendSize'
                LegSize = varargin{n+1};
                otherwise
                % We should probably just say that varargin was not a
                % recongnized option.... or this might get crazy...
                error([varargin{n},' was not found.'])
        end
    end
    HMarg = Marg(1);
    VMarg = Marg(2);
    LMarg = Marg(3);
    LWh = Wh(1);
    RWh = Wh(2);
    BWh = Wh(3);
    TWh = Wh(4);
    HLeg = LegSize(1);
    VLeg = LegSize(2);
    if ~isempty(LegPos)
        if strcmp(LegPos,'southeast')
            YLeg = BWh;
        elseif strcmp(LegPos,'east')
            YLeg = 0.5 - VLeg/2;
        elseif strcmp(LegPos,'northeast')
            YLeg = 1 - TWh - VLeg;
        end
    end
    pos = zeros(sum(Row),4);
    % 1 = TopWhiteSpace + length(Row)*Height + (length(Row)-1)*VMargin + BottomWhiteSpace
    % ==> H = (1 - (TopWhiteSpace + (length(Row)-1)*VMargin + BottomWhiteSpace)/ length(Row) 
    if isempty(LegPos)
    H = (1 - (TWh + (length(Row)-1)*VMarg + BWh))/ length(Row);
    for r = 1:length(Row)
        % 1 = LeftWhiteSpace + Row*Width + (Row-1)*HMargin + RightWhiteSpace
        %==> W = (1 - LeftWhiteSpace - (Row-1)*HMargin - RightWhiteSpace)/ Row
        W = (1 - LWh - (Row(r)-1)*HMarg - RWh)/ Row(r);
        for i = 1: Row(r)
            pos(i+sum(Row(1:r-1)),1) = LWh + (i-1)* W + (i-1)*HMarg;
            pos(i+sum(Row(1:r-1)),2) = BWh + (length(Row)-r)* H + (length(Row)-r)*VMarg;
            pos(i+sum(Row(1:r-1)),3) = W;
            pos(i+sum(Row(1:r-1)),4) = H;
        end
    end
    Leg = [];
    else
        H = (1 - (TWh + (length(Row)-1)*VMarg + BWh))/ length(Row);
    for r = 1:length(Row)
        % 1 = LeftWhiteSpace + Row*Width + (Row-1)*HMargin + LegendWidth + RightWhiteSpace
        %==> W = (1 - LeftWhiteSpace - (Row-1)*HMargin - LegendWidth - RightWhiteSpace)/ Row
        W = (1 - LWh - (Row(r)-1)*HMarg - HLeg - RWh - LMarg)/ Row(r);
        for i = 1: Row(r)
            pos(i+sum(Row(1:r-1)),1) = LWh + (i-1)* W + (i-1)*HMarg;
            pos(i+sum(Row(1:r-1)),2) = BWh + (length(Row)-r)* H + (length(Row)-r)*VMarg;
            pos(i+sum(Row(1:r-1)),3) = W;
            pos(i+sum(Row(1:r-1)),4) = H;
        end
    end
    XLeg = pos(end,1) + W + LMarg;
    Leg = [XLeg YLeg HLeg VLeg];
    end
    varargout{1} = pos;
    if nargout>1
        varargout{2} = Leg;
    end
        
end

