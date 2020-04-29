function Sol = rkCoefNonLinearRelations(xtry)

    Iterations = 1e4;
    vars = 36;
    if exist('xtry','var')
        xguess = xtry;
    else
        xguess = -1 + 2*rand(vars,1);
    end
    
    stages = 8;
    
    chP = chebyPoly(stages);
    Aeq = zeros(stages,vars);
    gg = 1;
    for m = 1:stages
        Aeq(m,gg:gg+m-1) = 1;
        gg = gg + m;
    end
    beq = [chP(2:end);1];
    optsCon = optimoptions('fmincon','MaxFunEvals',Iterations,...
        'MaxIter',Iterations,'TolFun',eps(0),'TolX',eps(0),...
        'TolCon',eps(.25));
    
    g = waitbar(0,'Solving function');
    setGlobalNumStep(1);
    [xs] = fmincon(@(x)funEval(x,Iterations),xguess,[],[],Aeq,beq,[],[],[],optsCon);
    
%     opts = optimoptions('fsolve','TolFun',eps,'TolX',eps,...
%         'MaxFunEvals',Iterations,'MaxIter',Iterations);
%     [xs,feval] = fsolve(@(x)funEval(x,Iterations),xguess,opts);
    
    waitbar(1);
    [feval,A,b,c] = funEval(xs);
    Sol.A = A;
    Sol.b = b;
    Sol.c = c;
    Sol.feval = feval;
    Sol.xs = xs;
    delete(g);

end

function setGlobalNumStep(val)
    global numSteps
    numSteps = val;
end
function r = getnumSteps
    global numSteps
    r = numSteps;
end
function chP = chebyPoly(stages)
    if stages == 1
        chP = 0;
    elseif stages == 2
        chP = [0;1];
    else
        chP = zeros(stages,1);
        chP(1) = 0;
        chP(end) = 1;
        chP(2:end-1) = cos((2.*((1*(stages-2)):-1:1) - 1)./(2*(stages-2)) * pi)/2 + 1/2;
    end
end

function [eval,A,b,c] = funEval(x,Iterations)

    if exist('Iterations','var')
        inSolver = 1;
        waitbar(getnumSteps/Iterations);
        setGlobalNumStep(getnumSteps+1);
    else
        inSolver = 0;
    end
    order = 6;
    
%     a0201 = x(1);
%     a0301 = x(2);
%     a0302 = x(3);
%     a0401 = x(4);
%     a0402 = x(5);
%     a0403 = x(6);
%     a0501 = x(7);
%     a0502 = x(8);
%     a0503 = x(9);
%     a0504 = x(10);
%     a0601 = x(11);
%     a0602 = x(12);
%     a0603 = x(13);
%     a0604 = x(14);
%     a0605 = x(15);
%     a0701 = x(16);
%     a0702 = x(17);
%     a0703 = x(18);
%     a0704 = x(19);
%     a0705 = x(20);
%     a0706 = x(21);
%     a0801 = x(22);
%     a0802 = x(23);
%     a0803 = x(24);
%     a0804 = x(25);
%     a0805 = x(26);
%     a0806 = x(27);
%     a0807 = x(28);
%     a0901 = x(29);
%     a0902 = x(30);
%     a0903 = x(31);
%     a0904 = x(32);
%     a0905 = x(33);
%     a0906 = x(34);
%     a0907 = x(35);
%     a0908 = x(36);
%     a1001 = x(37);
%     a1002 = x(38);
%     a1003 = x(39);
%     a1004 = x(40);
%     a1005 = x(41);
%     a1006 = x(42);
%     a1007 = x(43);
%     a1008 = x(44);
%     a1009 = x(45);
%     a1101 = x(46);
%     a1102 = x(47);
%     a1103 = x(48);
%     a1104 = x(49);
%     a1105 = x(50);
%     a1106 = x(51);
%     a1107 = x(52);
%     a1108 = x(53);
%     a1109 = x(54);
%     a1110 = x(55);
%     a1201 = x(56);
%     a1202 = x(57);
%     a1203 = x(58);
%     a1204 = x(59);
%     a1205 = x(60);
%     a1206 = x(61);
%     a1207 = x(62);
%     a1208 = x(63);
%     a1209 = x(64);
%     a1210 = x(65);
%     a1211 = x(66);
%     a1301 = x(67);
%     a1302 = x(68);
%     a1303 = x(69);
%     a1304 = x(70);
%     a1305 = x(71);
%     a1306 = x(72);
%     a1307 = x(73);
%     a1308 = x(74);
%     a1309 = x(75);
%     a1310 = x(76);
%     a1311 = x(77);
%     a1312 = x(78);
%     a1401 = x(79);
%     a1402 = x(80);
%     a1403 = x(81);
%     a1404 = x(82);
%     a1405 = x(83);
%     a1406 = x(84);
%     a1407 = x(85);
%     a1408 = x(86);
%     a1409 = x(87);
%     a1410 = x(88);
%     a1411 = x(89);
%     a1412 = x(90);
%     
%     b1 = x(91); b2 = x(92); b3 = x(93);
%     b4 = x(94); b5 = x(95); b6 = x(96);
%     b7 = x(97); b8 = x(98); b9 = x(99);
%     b10 = x(100); b11 = x(101); b12 = x(102); 
%     b13 = x(103); b14 = x(104);
% 
%     A = [
%         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%         a0201 0     0     0     0     0     0     0     0     0     0     0     0     0     0
%         a0301 a0302 0     0     0     0     0     0     0     0     0     0     0     0     0
%         a0401 a0402 a0403 0     0     0     0     0     0     0     0     0     0     0     0
%         a0501 a0502 a0503 a0504 0     0     0     0     0     0     0     0     0     0     0
%         a0601 a0602 a0603 a0604 a0605 0     0     0     0     0     0     0     0     0     0
%         a0701 a0702 a0703 a0704 a0705 a0706 0     0     0     0     0     0     0     0     0
%         a0801 a0802 a0803 a0804 a0805 a0806 a0807 0     0     0     0     0     0     0     0
%         a0901 a0902 a0903 a0904 a0905 a0906 a0907 a0908 0     0     0     0     0     0     0
%         a1001 a1002 a1003 a1004 a1005 a1006 a1007 a1008 a1009 0     0     0     0     0     0
%         a1101 a1102 a1103 a1104 a1105 a1106 a1107 a1108 a1109 a1110 0     0     0     0     0
%         a1201 a1202 a1203 a1204 a1205 a1206 a1207 a1208 a1209 a1210 a1211 0     0     0     0
%         a1301 a1302 a1303 a1304 a1305 a1306 a1307 a1308 a1309 a1310 a1311 a1312 0     0     0
%         a1401 a1402 a1403 a1404 a1405 a1406 a1407 a1408 a1409 a1410 a1411 a1412 a1412 0     0
%         b1    b2    b3    b4    b5    b6    b7    b8    b9    b10   b11   b12   b13   b14   0
%         ];
%     b1 = [
%         b1    b2    b3    b4    b5    b6    b7    b8    b9    b10   b11   b12   b13   b14   0
%         ];
%     b2 = x(105:119).';

    a21 = x(1);
    a31 = x(2);
    a32 = x(3);
    a41 = x(4);
    a42 = x(5);
    a43 = x(6);
    a51 = x(7);
    a52 = x(8);
    a53 = x(9);
    a54 = x(10);
    a61 = x(11);
    a62 = x(12);
    a63 = x(13);
    a64 = x(14);
    a65 = x(15);
    a71 = x(16);
    a72 = x(17);
    a73 = x(18);
    a74 = x(19);
    a75 = x(20);
    a76 = x(21);
    
    b1=x(22);b2=x(23);b3=x(24);b4=x(25);b5=x(26);b6=x(27);b7=x(28);
    A = [
        0   0   0   0   0   0   0   0
        a21 0   0   0   0   0   0   0
        a31 a32 0   0   0   0   0   0
        a41 a42 a43 0   0   0   0   0
        a51 a52 a53 a54 0   0   0   0
        a61 a62 a63 a64 a65 0   0   0
        a71 a72 a73 a74 a75 a76 0   0
        b1  b2  b3  b4  b5  b6  b7  0
        ];
    b1 = [
        b1  b2  b3  b4  b5  b6  b7  0];
    b2 = [x(29:36).'];
    c = sum(A,2);
    
    eval1 = oc_butcher(A,b1,c,order);
    eval2 = oc_butcher(A,b2,c,order-1);
    eval = [eval1;eval2];
    b = [b1;b2];
    if inSolver
        eval = norm(eval,1);
    end

end