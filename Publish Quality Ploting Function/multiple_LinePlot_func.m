function [ ] = multiple_LinePlot_func( varargin )
% This Function Plots multiple lines based on Elsevior journal standards.
    FigName='1';
    PlotCat=1;
    Size = 1;
    MarFer = 1;
    Marg = [0.05 0.05 0.01];
    Wh = [0.15 0.05 0.15 0.1];
    FigHei = 10;
    RowCount = 1;
    %     Subplot.XLabel = '';
    %     Subplot.YLabel = '';
    %     subplot.Legend = '';
    %     Subplot.Title = '';
    Subplot.LegLoc = 'NorthWest';
    Title = '';
    LegLoc = '';
    Leg = '';

    Subplot = varargin{1};
    for n = 2:2:numel(varargin)
        switch varargin{n}
            case 'FigureName'
                FigName = varargin{n+1};
            case 'Print'
                FigPrint = varargin{n+1};
            case 'Size'
                Size = varargin{n+1};
            case 'MarkerFrequency'
                MarFer = varargin{n+1};
            case 'CommonLegend'
                Leg = varargin{n+1};
            case 'CommonLegendLocation'
                LegLoc = varargin{n+1};
            case 'CommonLegendSize'
                LegSize = varargin{n+1};
            case 'CommonTitle'
                Title = varargin{n+1};
            case 'PerRow'
                RowCount = varargin{n+1};
            case 'FigureHeight'
                FigHei = varargin{n+1};
            case 'Margin'
                Marg = varargin{n+1};
            case 'WhiteSpace'
                Wh = varargin{n+1};
            otherwise
                % We should probably just say that varargin was not a
                % recongnized option.... or this might get crazy...
                error([varargin{n},' was not found.'])
        end
    end
    if Size == 2
        FigWid = 19;
    elseif Size == 1.5
        FigWid = 14;
    else
        FigWid = 9;
    end
    font_opt={'FontUnits','points','FontWeight','normal','FontSize',7,'FontName','Arial'};
    LaTeX_opt = {'interpreter','latex'};
    FigSize = {'Units', 'centimeters', 'Position', [-30, 3,FigWid, FigHei], 'Paperpositionmode','auto'};
    plot_Line_opt = {'LineWidth',1};
    shape1= {'-sk' '-ok' '-vk' '-^k' '-*k' '-xk' '-dk' '-pk' '-hk' '->k' '-<k' '--sk' '--ok' '--vk' '--^k' '--*k' '--xk' '--dk' '--pk' '--hk' '-->k' '--<k'};
    shape2= {'-sb' '-og' '-vr' '-^k' '-vc' '-ys' '-kd' '-r^' '-gv' '-b>' '-c<' '-mp' '-yh' '-sk' '-ok' '-vk' '-^k' '-*k' '-xk' '-dk' '-pk' '-hk' '->k' '-<k' '--sk' '--ok' '--vk' '--^k' '--*k' '--xk'};

    if strcmp(FigPrint,'Color')
        PlotShape=shape2;
        PlotCat = 2;
    elseif strcmp(FigPrint,'Black')
        PlotShape=shape1;
    else
        PlotShape = FigPrint;
        PlotCat =2;
    end
    PlotShape1 = PlotShape;
    for i = 1:length(PlotShape)
        PlotShape{i}(regexp(PlotShape{i},'[F]'))=[];
    end


    figure('Name',FigName) %Steel
    set(gcf, FigSize{:});
    p = uipanel('Parent',gcf,'BorderType','none');
    set(p,'Title', Title, 'TitlePosition', 'centertop',font_opt{:});
    clf
    movegui(gcf,'center');
    [pos LegPos] = subplot_position_manager('PerRow',RowCount,'Margin',Marg,'WhiteSpace',Wh,'LegendPosition',LegLoc);
    row = 1;
    for kk = 1:sum(RowCount)
        legend_opt= {'Location',Subplot(kk).LegLoc};
        if (kk > sum(RowCount(1:row))),  row = row + 1; end
        subplotloc = max(RowCount)*(row-1) + (kk - sum(RowCount(1:row-1))-1)*round(max(RowCount)/RowCount(row))+1;
        handle(kk,:) = axes; % subplot(length(RowCount),max(RowCount),subplotloc);
        for m = 1:size(Subplot(kk).Y,2)
            h(m,:) = plot(Subplot(kk).X(:,m),Subplot(kk).Y(:,m),PlotShape{m});
            grid on
            set(h(m,:),plot_Line_opt{:});
            if PlotCat == 1
                h(m,:).MarkerSize = 5;
                h(m,:).MarkerFaceColor = 'k';
            end
            if MarFer ~= 1;
                h(m,:).MarkerIndices = m:(MarFer):size(Subplot(kk).Y,1);
            end
            if any(~cellfun('isempty',strfind(PlotShape1(m),'F')))
                set(h(m,:),'MarkerFaceColor','k');
            end
            axis(Subplot(kk).Axis(1:4))
            hold on
        end
        if ~isempty(Subplot(kk).XLabel), xlabel(Subplot(kk).XLabel,font_opt{:},LaTeX_opt{:}); end
        if length(Subplot(kk).Axis) >4
            if (Subplot(kk).Axis(7)==0), set(handle(kk,:),'XTickLabel',[]);end
        end
        if ~isempty(Subplot(kk).YLabel), ylabel(Subplot(kk).YLabel,font_opt{:},LaTeX_opt{:}); end
        if length(Subplot(kk).Axis) >4
            if (Subplot(kk).Axis(8)==0), set(handle(kk,:),'YTickLabel',[]);end
        end
        if ~isempty(Subplot(kk).Legend), legend(Subplot(kk).Legend,legend_opt{:},font_opt{:},LaTeX_opt{:}); end
        if ~isempty(Subplot(kk).Title), title(Subplot(kk).Title,font_opt{:},LaTeX_opt{:}); end
    end
    if ~strcmp(Leg,'')
        hL = legend(Leg,font_opt{:},LaTeX_opt{:});
        CurLegLoc = get(hL,'Position');
        if (CurLegLoc(3)> LegPos(3)*1.2 | CurLegLoc(3)< LegPos(3)*0.8 | CurLegLoc(4)> LegPos(4)*1.2 | CurLegLoc(4)< LegPos(4)*0.8)
            LegSize = CurLegLoc(3:4);
            [pos LegPos] = subplot_position_manager('PerRow',RowCount,'Margin',Marg,'WhiteSpace',Wh,'LegendPosition',LegLoc,'LegendSize',LegSize);
        end
        set(hL,'Position', LegPos,'Units', 'normalized');
    end
    
    for kk = 1:sum(RowCount)
        if (Subplot(kk).Axis(2)-Subplot(kk).Axis(1))< 1
           XTICK = round(Subplot(kk).Axis(1):(Subplot(kk).Axis(2)-Subplot(kk).Axis(1))/Subplot(kk).Axis(5):Subplot(kk).Axis(2),3);
        elseif (Subplot(kk).Axis(2)-Subplot(kk).Axis(1))< 10
            XTICK = round(Subplot(kk).Axis(1):(Subplot(kk).Axis(2)-Subplot(kk).Axis(1))/Subplot(kk).Axis(5):Subplot(kk).Axis(2),2);
        elseif (Subplot(kk).Axis(2)-Subplot(kk).Axis(1))< 100
            XTICK = round(Subplot(kk).Axis(1):(Subplot(kk).Axis(2)-Subplot(kk).Axis(1))/Subplot(kk).Axis(5):Subplot(kk).Axis(2),1);
        else
            XTICK = round(Subplot(kk).Axis(1):(Subplot(kk).Axis(2)-Subplot(kk).Axis(1))/Subplot(kk).Axis(5):Subplot(kk).Axis(2));
        end
        
        if (Subplot(kk).Axis(4)-Subplot(kk).Axis(3))< 1
           YTICK = round(Subplot(kk).Axis(3):(Subplot(kk).Axis(4)-Subplot(kk).Axis(3))/Subplot(kk).Axis(6):Subplot(kk).Axis(4),3);
        elseif (Subplot(kk).Axis(4)-Subplot(kk).Axis(3))< 10
            YTICK = round(Subplot(kk).Axis(3):(Subplot(kk).Axis(4)-Subplot(kk).Axis(3))/Subplot(kk).Axis(6):Subplot(kk).Axis(4),2);
        elseif (Subplot(kk).Axis(4)-Subplot(kk).Axis(3))< 100
            YTICK = round(Subplot(kk).Axis(3):(Subplot(kk).Axis(4)-Subplot(kk).Axis(3))/Subplot(kk).Axis(6):Subplot(kk).Axis(4),1);
        else
            YTICK = round(Subplot(kk).Axis(3):(Subplot(kk).Axis(4)-Subplot(kk).Axis(3))/Subplot(kk).Axis(6):Subplot(kk).Axis(4));
        end
            
        set(handle(kk,:),'Units','normalized','XTick',XTICK,'YTick',YTICK,font_opt{:});%,'yticklabel',num2str(get(handle(kk,:),'ytick')','%g'));
%         ytickformat('%3.2f')
        %         xtickformat('%.2f')
        set(handle(kk,:),'position',pos(kk,:))
    end
%     if ~strcmp(Leg,'')
%         hL = legend(Leg,font_opt{:},LaTeX_opt{:});
%         CurLegLoc = get(hL,'Position');
%         if (CurLegLoc(3)> LegPos(3)*1.2 | CurLegLoc(3)< LegPos(3)*0.8 | CurLegLoc(4)> LegPos(4)*1.2 | CurLegLoc(4)< LegPos(4)*0.8)
%             LegSize = CurLegLoc(3:4);
%             [pos LegPos] = subplot_position_manager('PerRow',RowCount,'Margin',Marg,'WhiteSpace',Wh,'LegendPosition',LegLoc,'LegendSize',LegSize);
%             row = 1;
%             clf
%             for kk = 1:sum(RowCount)
%                 legend_opt= {'Location',Subplot(kk).LegLoc};
%                 if (kk > sum(RowCount(1:row))),  row = row + 1; end
%                 subplotloc = max(RowCount)*(row-1) + (kk - sum(RowCount(1:row-1))-1)*(max(RowCount)/RowCount(row))+1;
%                 posloc = max(RowCount)*(row-1) + (kk - sum(RowCount(1:row-1)));
%                 handle(kk,:) = subplot(length(RowCount),max(RowCount),subplotloc);
%                 for m = 1:size(Subplot(kk).Y,2)
%                     h(m,:) = plot(Subplot(kk).X(:,m),Subplot(kk).Y(:,m),PlotShape{m});
%                     grid on
%                     set(h(m,:),plot_Line_opt{:});
%                     if PlotCat == 1
%                         h(m,:).MarkerSize = 5;
%                         h(m,:).MarkerFaceColor = 'k';
%                     end
%                     if MarFer ~= 1;
%                         h(m,:).MarkerIndices = m:(MarFer):size(Subplot(kk).Y,1);
%                     end
%                     if any(~cellfun('isempty',strfind(PlotShape1(m),'F')))
%                         set(h(m,:),'MarkerFaceColor','k');
%                     end
%                     axis(Subplot(kk).Axis(1:4))
%                     hold on
%                 end
%             end
%             for kk = 1:sum(RowCount)
%                 set(handle(kk,:),'Units','normalized','XTick',Subplot(kk).Axis(1):(Subplot(kk).Axis(2)-Subplot(kk).Axis(1))/Subplot(kk).Axis(5):Subplot(kk).Axis(2),'YTick',Subplot(kk).Axis(3):(Subplot(kk).Axis(4)-Subplot(kk).Axis(3))/Subplot(kk).Axis(6):Subplot(kk).Axis(4),font_opt{:});
%                 ytickformat('%.2f')
%                 %         xtickformat('%.2f')
%                 set(handle(kk,:),'position',pos(posloc,:))
%                 if ~isempty(Subplot(kk).XLabel), xlabel(Subplot(kk).XLabel,font_opt{:},LaTeX_opt{:}); end
%                 if length(Subplot(kk).Axis) >4
%                     if (Subplot(kk).Axis(7)==0), set(handle(kk,:),'XTickLabel',[]);end
%                 end
%                 if ~isempty(Subplot(kk).YLabel), ylabel(Subplot(kk).YLabel,font_opt{:},LaTeX_opt{:}); end
%                 if length(Subplot(kk).Axis) >4
%                     if (Subplot(kk).Axis(8)==0), set(handle(kk,:),'YTickLabel',[]);end
%                 end
%                 if ~isempty(Subplot(kk).Legend), legend(Subplot(kk).Legend,legend_opt{:},font_opt{:},LaTeX_opt{:}); end
%                 if ~isempty(Subplot(kk).Title), title(Subplot(kk).Title,font_opt{:},LaTeX_opt{:}); end
%             end
%             hL = legend(Leg,font_opt{:},LaTeX_opt{:});
%         end
% 
%         set(hL,'Position', LegPos,'Units', 'normalized');
%     end
% 

    % Ask user if they wish to save the figure
    promptMessage = sprintf(['Do you wish to save this figure: ',Title,'?']);
    titleBarCaption = 'Yes or No';
    button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
    if strcmpi(button, 'Yes')
        buttonSelections = 1;
    else
        buttonSelections = 0;
    end
    if buttonSelections == 1
        [FileName,PathName,FilterIndex]= uiputfile(...
            {'*.eps','Vector Figures(*.eps)';...
            '*.pdf','PDF File(*.pdf)';...
            '*.fig','Matlab Figures (*.fig)';...
            '*.jpg;*.tif;*.png;*.gif','Other Image Files(*.jpg)'},'Save Image',[FigName]);
        if FilterIndex == 1
            print(figure(1),'-depsc2',[PathName,FileName]);
        elseif FilterIndex == 2
            print(figure(1),'-dpdf',[PathName,FileName]);
        elseif FilterIndex == 4
            print(figure(1),'-djpeg',[PathName,FileName]);
        elseif FilterIndex == 3
            savefig([PathName,FileName]);
        end
    end
end

