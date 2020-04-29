% What someone on Minecraft Owes me...
% Aug 24 start date
function WhatIsOwed()
    startDate = [2015 8 25 0 0 0];
    x0 = 64+12; % blocks...
    rate = 0.0001;

    days_on_Loan = etime(clock, startDate)/3600/24;
    MineCraft_Days_in_Hour = 3;
    MineCraft_Days_in_Day = MineCraft_Days_in_Hour*24;
    tend = days_on_Loan * MineCraft_Days_in_Day;

    owed = x0*exp(rate*tend);
    blocks = floor(owed);
    remainder = owed-blocks;
    bars = floor(remainder*9);
    remainder = 9*remainder - bars;
    nuggets = floor(remainder*9);
    remainder = 9*remainder-nuggets;

    fprintf([
            '\n Ammount Owed...\n'...
            '    Blocks = %i\n'...
            '      Bars = %i\n'... 
            '   Nuggets = %i\n'...
            ' Remainder = %0.3f\n'],...
            blocks,bars,nuggets,remainder); 
end