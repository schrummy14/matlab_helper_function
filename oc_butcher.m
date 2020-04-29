function coneq=oc_butcher(A,b,c,p)
    % function coneq=oc_butcher(A,b,c,p)
    %
    % Order conditions for RKMs.
    % This version is based on Butcher's approach.
    %
    % Assumes `p>1`.

    if isempty(A) 
        ordernumFlag = 1;
        A = 0;
        b = 0;
        c = 0;
    else
        ordernumFlag = 0;
    end

    if size(b,2)>1,b = b.';end
    if size(c,2)>1,c = c.';end

    coneq = [];
    if p>=2
        coneq(1,1)=c'*b-1/2;
    end

    if p>=3   
      coneq(2,1)=c'.^2*b-1/3;
      coneq(3,1)=b'*A*c-1/6;
    end

    if p>=4  
      coneq(4,1)=c'.^3*b-1/4;
      coneq(5,1)=(b'.*c')*A*c-1/8;
      coneq(6,1)=b'*A*c.^2-1/12;
      coneq(7,1)=b'*A^2*c-1/24;
    end

    if p>=5 
      coneq(8,1)=c'.^4*b-1/5;
      coneq(9,1)=(b.*c.^2)'*A*c-1/10;
      coneq(10,1)=b'*(A*c).^2-1/20;
      coneq(11,1)=(b.*c)'*A*c.^2-1/15;
      coneq(12,1)=b'*A*c.^3-1/20;
      coneq(13,1)=(b.*c)'*A^2*c-1/30;
      coneq(14,1)=b'*A*diag(c)*A*c-1/40;
      coneq(15,1)=b'*A^2*c.^2-1/60;
      coneq(16,1)=b'*A^3*c-1/120;
    end

    if p>=6
      coneq(17,1)=c'.^5*b-1/6;
      coneq(18,1)=b'*diag(c).^3*A*c-1/12;
      coneq(19,1)=b'*diag(c)*(A*c).^2-1/24;
      coneq(20,1)=b'*diag(c).^2*A*c.^2-1/18;
      coneq(21,1)=b'*((A*c.^2).*(A*c))-1/36;
      coneq(22,1)=b'*diag(c)*A*c.^3-1/24;
      coneq(23,1)=b'*A*c.^4-1/30;
      coneq(24,1)=b'*diag(c).^2*A^2*c-1/36;
      coneq(25,1)=b'*((A^2*c).*(A*c))-1/72;
      coneq(26,1)=b'*diag(c)*A*diag(c)*A*c-1/48;
      coneq(27,1)=b'*A*diag(c).^2*A*c-1/60;
      coneq(28,1)=b'*A*(A*c).^2-1/120;
      coneq(29,1)=b'*diag(c)*A^2*c.^2-1/72;
      coneq(30,1)=b'*A*diag(c)*A*c.^2-1/90;
      coneq(31,1)=b'*A^2*c.^3-1/120;
      coneq(32,1)=b'*diag(c)*A^3*c-1/144;
      coneq(33,1)=b'*A*diag(c)*A^2*c-1/180;
      coneq(34,1)=b'*A^2*diag(c)*A*c-1/240;
      coneq(35,1)=b'*A^3*c.^2-1/360;
      coneq(36,1)=b'*A^4*c-1/720;
    end

    if p>=7
      % order 7 conditions:
      coneq(37)=b'*((A*c.^5))-1/42;
      coneq(38)=b'*((A*(A*c.^4)))-1/210;
      coneq(39)=b'*((A*(A*(A*c.^3))))-1/840;
      coneq(40)=b'*((A*(A*(A*(A*c.^2)))))-1/2520;
      coneq(41)=b'*((A*(A*(A*(A*(A*c))))))-1/5040;
      coneq(42)=b'*((A*(A*(A*c.*(A*c)))))-1/1680;
      coneq(43)=b'*((A*(A*c.*(A*c.^2))))-1/630;
      coneq(44)=b'*((A*(A*c.*(A*(A*c)))))-1/1260;
      coneq(45)=b'*((A*(A*c.^2.*(A*c))))-1/420;
      coneq(46)=b'*((A*(A*(A*c).*(A*c))))-1/840;
      coneq(47)=b'*((A*c.*(A*c.^3)))-1/168;
      coneq(48)=b'*((A*c.*(A*(A*c.^2))))-1/504;
      coneq(49)=b'*((A*c.*(A*(A*(A*c)))))-1/1008;
      coneq(50)=b'*((A*c.*(A*c.*(A*c))))-1/336;
      coneq(51)=b'*((A*c.^2.*(A*c.^2)))-1/126;
      coneq(52)=b'*((A*c.^2.*(A*(A*c))))-1/252;
      coneq(53)=b'*((A*c.^3.*(A*c)))-1/84;
      coneq(54)=b'*((A*(A*c).*(A*c.^2)))-1/252;
      coneq(55)=b'*((A*(A*c).*(A*(A*c))))-1/504;
      coneq(56)=b'*((A*c.*(A*c).*(A*c)))-1/168;
      coneq(57)=b'*((A*c.^4).*c)-1/35;
      coneq(58)=b'*((A*(A*c.^3)).*c)-1/140;
      coneq(59)=b'*((A*(A*(A*c.^2))).*c)-1/420;
      coneq(60)=b'*((A*(A*(A*(A*c)))).*c)-1/840;
      coneq(61)=b'*((A*(A*c.*(A*c))).*c)-1/280;
      coneq(62)=b'*((A*c.*(A*c.^2)).*c)-1/105;
      coneq(63)=b'*((A*c.*(A*(A*c))).*c)-1/210;
      coneq(64)=b'*((A*c.^2.*(A*c)).*c)-1/70;
      coneq(65)=b'*((A*(A*c).*(A*c)).*c)-1/140;
      coneq(66)=b'*((A*c.^3).*c.^2)-1/28;
      coneq(67)=b'*((A*(A*c.^2)).*c.^2)-1/84;
      coneq(68)=b'*((A*(A*(A*c))).*c.^2)-1/168;
      coneq(69)=b'*((A*c.*(A*c)).*c.^2)-1/56;
      coneq(70)=b'*((A*c.^2).*c.^3)-1/21;
      coneq(71)=b'*((A*(A*c)).*c.^3)-1/42;
      coneq(72)=b'*((A*c).*c.^4)-1/14;
      coneq(73)=b'*((A*c.^2).*(A*c.^2))-1/63;
      coneq(74)=b'*((A*(A*c)).*(A*c.^2))-1/126;
      coneq(75)=b'*((A*(A*c)).*(A*(A*c)))-1/252;
      coneq(76)=b'*((A*c.^3).*(A*c))-1/56;
      coneq(77)=b'*((A*(A*c.^2)).*(A*c))-1/168;
      coneq(78)=b'*((A*(A*(A*c))).*(A*c))-1/336;
      coneq(79)=b'*((A*c.*(A*c)).*(A*c))-1/112;
      coneq(80)=b'*((A*c.^2).*(A*c).*c)-1/42;
      coneq(81)=b'*((A*(A*c)).*(A*c).*c)-1/84;
      coneq(82)=b'*((A*c).*(A*c).*c.^2)-1/28;
      coneq(83)=b'*((A*c).*(A*c).*(A*c))-1/56;
      coneq(84)=b'*(c.^6)-1/7;
    end

    if p>=8
      % order 8 conditions:
      coneq(85)=b'*((A*c.^6))-1/56;
      coneq(86)=b'*((A*(A*c.^5)))-1/336;
      coneq(87)=b'*((A*(A*(A*c.^4))))-1/1680;
      coneq(88)=b'*((A*(A*(A*(A*c.^3)))))-1/6720;
      coneq(89)=b'*((A*(A*(A*(A*(A*c.^2))))))-1/20160;
      coneq(90)=b'*((A*(A*(A*(A*(A*(A*c)))))))-1/40320;
      coneq(91)=b'*((A*(A*(A*(A*c.*(A*c))))))-1/13440;
      coneq(92)=b'*((A*(A*(A*c.*(A*c.^2)))))-1/5040;
      coneq(93)=b'*((A*(A*(A*c.*(A*(A*c))))))-1/10080;
      coneq(94)=b'*((A*(A*(A*c.^2.*(A*c)))))-1/3360;
      coneq(95)=b'*((A*(A*(A*(A*c).*(A*c)))))-1/6720;
      coneq(96)=b'*((A*(A*c.*(A*c.^3))))-1/1344;
      coneq(97)=b'*((A*(A*c.*(A*(A*c.^2)))))-1/4032;
      coneq(98)=b'*((A*(A*c.*(A*(A*(A*c))))))-1/8064;
      coneq(99)=b'*((A*(A*c.*(A*c.*(A*c)))))-1/2688;
      coneq(100)=b'*((A*(A*c.^2.*(A*c.^2))))-1/1008;
      coneq(101)=b'*((A*(A*c.^2.*(A*(A*c)))))-1/2016;
      coneq(102)=b'*((A*(A*c.^3.*(A*c))))-1/672;
      coneq(103)=b'*((A*(A*(A*c).*(A*c.^2))))-1/2016;
      coneq(104)=b'*((A*(A*(A*c).*(A*(A*c)))))-1/4032;
      coneq(105)=b'*((A*(A*c.*(A*c).*(A*c))))-1/1344;
      coneq(106)=b'*((A*c.*(A*c.^4)))-1/280;
      coneq(107)=b'*((A*c.*(A*(A*c.^3))))-1/1120;
      coneq(108)=b'*((A*c.*(A*(A*(A*c.^2)))))-1/3360;
      coneq(109)=b'*((A*c.*(A*(A*(A*(A*c))))))-1/6720;
      coneq(110)=b'*((A*c.*(A*(A*c.*(A*c)))))-1/2240;
      coneq(111)=b'*((A*c.*(A*c.*(A*c.^2))))-1/840;
      coneq(112)=b'*((A*c.*(A*c.*(A*(A*c)))))-1/1680;
      coneq(113)=b'*((A*c.*(A*c.^2.*(A*c))))-1/560;
      coneq(114)=b'*((A*c.*(A*(A*c).*(A*c))))-1/1120;
      coneq(115)=b'*((A*c.^2.*(A*c.^3)))-1/224;
      coneq(116)=b'*((A*c.^2.*(A*(A*c.^2))))-1/672;
      coneq(117)=b'*((A*c.^2.*(A*(A*(A*c)))))-1/1344;
      coneq(118)=b'*((A*c.^2.*(A*c.*(A*c))))-1/448;
      coneq(119)=b'*((A*c.^3.*(A*c.^2)))-1/168;
      coneq(120)=b'*((A*c.^3.*(A*(A*c))))-1/336;
      coneq(121)=b'*((A*c.^4.*(A*c)))-1/112;
      coneq(122)=b'*((A*(A*c.^2).*(A*c.^2)))-1/504;
      coneq(123)=b'*((A*(A*c.^2).*(A*(A*c))))-1/1008;
      coneq(124)=b'*((A*(A*(A*c)).*(A*(A*c))))-1/2016;
      coneq(125)=b'*((A*(A*c).*(A*c.^3)))-1/448;
      coneq(126)=b'*((A*(A*c).*(A*(A*c.^2))))-1/1344;
      coneq(127)=b'*((A*(A*c).*(A*(A*(A*c)))))-1/2688;
      coneq(128)=b'*((A*(A*c).*(A*c.*(A*c))))-1/896;
      coneq(129)=b'*((A*c.*(A*c).*(A*c.^2)))-1/336;
      coneq(130)=b'*((A*c.*(A*c).*(A*(A*c))))-1/672;
      coneq(131)=b'*((A*c.^2.*(A*c).*(A*c)))-1/224;
      coneq(132)=b'*((A*(A*c).*(A*c).*(A*c)))-1/448;
      coneq(133)=b'*((A*c.^5).*c)-1/48;
      coneq(134)=b'*((A*(A*c.^4)).*c)-1/240;
      coneq(135)=b'*((A*(A*(A*c.^3))).*c)-1/960;
      coneq(136)=b'*((A*(A*(A*(A*c.^2)))).*c)-1/2880;
      coneq(137)=b'*((A*(A*(A*(A*(A*c))))).*c)-1/5760;
      coneq(138)=b'*((A*(A*(A*c.*(A*c)))).*c)-1/1920;
      coneq(139)=b'*((A*(A*c.*(A*c.^2))).*c)-1/720;
      coneq(140)=b'*((A*(A*c.*(A*(A*c)))).*c)-1/1440;
      coneq(141)=b'*((A*(A*c.^2.*(A*c))).*c)-1/480;
      coneq(142)=b'*((A*(A*(A*c).*(A*c))).*c)-1/960;
      coneq(143)=b'*((A*c.*(A*c.^3)).*c)-1/192;
      coneq(144)=b'*((A*c.*(A*(A*c.^2))).*c)-1/576;
      coneq(145)=b'*((A*c.*(A*(A*(A*c)))).*c)-1/1152;
      coneq(146)=b'*((A*c.*(A*c.*(A*c))).*c)-1/384;
      coneq(147)=b'*((A*c.^2.*(A*c.^2)).*c)-1/144;
      coneq(148)=b'*((A*c.^2.*(A*(A*c))).*c)-1/288;
      coneq(149)=b'*((A*c.^3.*(A*c)).*c)-1/96;
      coneq(150)=b'*((A*(A*c).*(A*c.^2)).*c)-1/288;
      coneq(151)=b'*((A*(A*c).*(A*(A*c))).*c)-1/576;
      coneq(152)=b'*((A*c.*(A*c).*(A*c)).*c)-1/192;
      coneq(153)=b'*((A*c.^4).*c.^2)-1/40;
      coneq(154)=b'*((A*(A*c.^3)).*c.^2)-1/160;
      coneq(155)=b'*((A*(A*(A*c.^2))).*c.^2)-1/480;
      coneq(156)=b'*((A*(A*(A*(A*c)))).*c.^2)-1/960;
      coneq(157)=b'*((A*(A*c.*(A*c))).*c.^2)-1/320;
      coneq(158)=b'*((A*c.*(A*c.^2)).*c.^2)-1/120;
      coneq(159)=b'*((A*c.*(A*(A*c))).*c.^2)-1/240;
      coneq(160)=b'*((A*c.^2.*(A*c)).*c.^2)-1/80;
      coneq(161)=b'*((A*(A*c).*(A*c)).*c.^2)-1/160;
      coneq(162)=b'*((A*c.^3).*c.^3)-1/32;
      coneq(163)=b'*((A*(A*c.^2)).*c.^3)-1/96;
      coneq(164)=b'*((A*(A*(A*c))).*c.^3)-1/192;
      coneq(165)=b'*((A*c.*(A*c)).*c.^3)-1/64;
      coneq(166)=b'*((A*c.^2).*c.^4)-1/24;
      coneq(167)=b'*((A*(A*c)).*c.^4)-1/48;
      coneq(168)=b'*((A*c).*c.^5)-1/16;
      coneq(169)=b'*((A*c.^3).*(A*c.^2))-1/96;
      coneq(170)=b'*((A*(A*c.^2)).*(A*c.^2))-1/288;
      coneq(171)=b'*((A*(A*(A*c))).*(A*c.^2))-1/576;
      coneq(172)=b'*((A*c.*(A*c)).*(A*c.^2))-1/192;
      coneq(173)=b'*((A*c.^3).*(A*(A*c)))-1/192;
      coneq(174)=b'*((A*(A*c.^2)).*(A*(A*c)))-1/576;
      coneq(175)=b'*((A*(A*(A*c))).*(A*(A*c)))-1/1152;
      coneq(176)=b'*((A*c.*(A*c)).*(A*(A*c)))-1/384;
      coneq(177)=b'*((A*c.^4).*(A*c))-1/80;
      coneq(178)=b'*((A*(A*c.^3)).*(A*c))-1/320;
      coneq(179)=b'*((A*(A*(A*c.^2))).*(A*c))-1/960;
      coneq(180)=b'*((A*(A*(A*(A*c)))).*(A*c))-1/1920;
      coneq(181)=b'*((A*(A*c.*(A*c))).*(A*c))-1/640;
      coneq(182)=b'*((A*c.*(A*c.^2)).*(A*c))-1/240;
      coneq(183)=b'*((A*c.*(A*(A*c))).*(A*c))-1/480;
      coneq(184)=b'*((A*c.^2.*(A*c)).*(A*c))-1/160;
      coneq(185)=b'*((A*(A*c).*(A*c)).*(A*c))-1/320;
      coneq(186)=b'*((A*c.^2).*(A*c.^2).*c)-1/72;
      coneq(187)=b'*((A*(A*c)).*(A*c.^2).*c)-1/144;
      coneq(188)=b'*((A*(A*c)).*(A*(A*c)).*c)-1/288;
      coneq(189)=b'*((A*c.^3).*(A*c).*c)-1/64;
      coneq(190)=b'*((A*(A*c.^2)).*(A*c).*c)-1/192;
      coneq(191)=b'*((A*(A*(A*c))).*(A*c).*c)-1/384;
      coneq(192)=b'*((A*c.*(A*c)).*(A*c).*c)-1/128;
      coneq(193)=b'*((A*c.^2).*(A*c).*c.^2)-1/48;
      coneq(194)=b'*((A*(A*c)).*(A*c).*c.^2)-1/96;
      coneq(195)=b'*((A*c).*(A*c).*c.^3)-1/32;
      coneq(196)=b'*((A*c.^2).*(A*c).*(A*c))-1/96;
      coneq(197)=b'*((A*(A*c)).*(A*c).*(A*c))-1/192;
      coneq(198)=b'*((A*c).*(A*c).*(A*c).*c)-1/64;
      coneq(199)=b'*(c.^7)-1/8;
    end

    if p>=9
      % order 9 conditions:
      coneq(200)=b'*((A*c.^7))-1/72;
      coneq(201)=b'*((A*(A*c.^6)))-1/504;
      coneq(202)=b'*((A*(A*(A*c.^5))))-1/3024;
      coneq(203)=b'*((A*(A*(A*(A*c.^4)))))-1/15120;
      coneq(204)=b'*((A*(A*(A*(A*(A*c.^3))))))-1/60480;
      coneq(205)=b'*((A*(A*(A*(A*(A*(A*c.^2)))))))-1/181440;
      coneq(206)=b'*((A*(A*(A*(A*(A*(A*(A*c))))))))-1/362880;
      coneq(207)=b'*((A*(A*(A*(A*(A*c.*(A*c)))))))-1/120960;
      coneq(208)=b'*((A*(A*(A*(A*c.*(A*c.^2))))))-1/45360;
      coneq(209)=b'*((A*(A*(A*(A*c.*(A*(A*c)))))))-1/90720;
      coneq(210)=b'*((A*(A*(A*(A*c.^2.*(A*c))))))-1/30240;
      coneq(211)=b'*((A*(A*(A*(A*(A*c).*(A*c))))))-1/60480;
      coneq(212)=b'*((A*(A*(A*c.*(A*c.^3)))))-1/12096;
      coneq(213)=b'*((A*(A*(A*c.*(A*(A*c.^2))))))-1/36288;
      coneq(214)=b'*((A*(A*(A*c.*(A*(A*(A*c)))))))-1/72576;
      coneq(215)=b'*((A*(A*(A*c.*(A*c.*(A*c))))))-1/24192;
      coneq(216)=b'*((A*(A*(A*c.^2.*(A*c.^2)))))-1/9072;
      coneq(217)=b'*((A*(A*(A*c.^2.*(A*(A*c))))))-1/18144;
      coneq(218)=b'*((A*(A*(A*c.^3.*(A*c)))))-1/6048;
      coneq(219)=b'*((A*(A*(A*(A*c).*(A*c.^2)))))-1/18144;
      coneq(220)=b'*((A*(A*(A*(A*c).*(A*(A*c))))))-1/36288;
      coneq(221)=b'*((A*(A*(A*c.*(A*c).*(A*c)))))-1/12096;
      coneq(222)=b'*((A*(A*c.*(A*c.^4))))-1/2520;
      coneq(223)=b'*((A*(A*c.*(A*(A*c.^3)))))-1/10080;
      coneq(224)=b'*((A*(A*c.*(A*(A*(A*c.^2))))))-1/30240;
      coneq(225)=b'*((A*(A*c.*(A*(A*(A*(A*c)))))))-1/60480;
      coneq(226)=b'*((A*(A*c.*(A*(A*c.*(A*c))))))-1/20160;
      coneq(227)=b'*((A*(A*c.*(A*c.*(A*c.^2)))))-1/7560;
      coneq(228)=b'*((A*(A*c.*(A*c.*(A*(A*c))))))-1/15120;
      coneq(229)=b'*((A*(A*c.*(A*c.^2.*(A*c)))))-1/5040;
      coneq(230)=b'*((A*(A*c.*(A*(A*c).*(A*c)))))-1/10080;
      coneq(231)=b'*((A*(A*c.^2.*(A*c.^3))))-1/2016;
      coneq(232)=b'*((A*(A*c.^2.*(A*(A*c.^2)))))-1/6048;
      coneq(233)=b'*((A*(A*c.^2.*(A*(A*(A*c))))))-1/12096;
      coneq(234)=b'*((A*(A*c.^2.*(A*c.*(A*c)))))-1/4032;
      coneq(235)=b'*((A*(A*c.^3.*(A*c.^2))))-1/1512;
      coneq(236)=b'*((A*(A*c.^3.*(A*(A*c)))))-1/3024;
      coneq(237)=b'*((A*(A*c.^4.*(A*c))))-1/1008;
      coneq(238)=b'*((A*(A*(A*c.^2).*(A*c.^2))))-1/4536;
      coneq(239)=b'*((A*(A*(A*c.^2).*(A*(A*c)))))-1/9072;
      coneq(240)=b'*((A*(A*(A*(A*c)).*(A*(A*c)))))-1/18144;
      coneq(241)=b'*((A*(A*(A*c).*(A*c.^3))))-1/4032;
      coneq(242)=b'*((A*(A*(A*c).*(A*(A*c.^2)))))-1/12096;
      coneq(243)=b'*((A*(A*(A*c).*(A*(A*(A*c))))))-1/24192;
      coneq(244)=b'*((A*(A*(A*c).*(A*c.*(A*c)))))-1/8064;
      coneq(245)=b'*((A*(A*c.*(A*c).*(A*c.^2))))-1/3024;
      coneq(246)=b'*((A*(A*c.*(A*c).*(A*(A*c)))))-1/6048;
      coneq(247)=b'*((A*(A*c.^2.*(A*c).*(A*c))))-1/2016;
      coneq(248)=b'*((A*(A*(A*c).*(A*c).*(A*c))))-1/4032;
      coneq(249)=b'*((A*c.*(A*c.^5)))-1/432;
      coneq(250)=b'*((A*c.*(A*(A*c.^4))))-1/2160;
      coneq(251)=b'*((A*c.*(A*(A*(A*c.^3)))))-1/8640;
      coneq(252)=b'*((A*c.*(A*(A*(A*(A*c.^2))))))-1/25920;
      coneq(253)=b'*((A*c.*(A*(A*(A*(A*(A*c)))))))-1/51840;
      coneq(254)=b'*((A*c.*(A*(A*(A*c.*(A*c))))))-1/17280;
      coneq(255)=b'*((A*c.*(A*(A*c.*(A*c.^2)))))-1/6480;
      coneq(256)=b'*((A*c.*(A*(A*c.*(A*(A*c))))))-1/12960;
      coneq(257)=b'*((A*c.*(A*(A*c.^2.*(A*c)))))-1/4320;
      coneq(258)=b'*((A*c.*(A*(A*(A*c).*(A*c)))))-1/8640;
      coneq(259)=b'*((A*c.*(A*c.*(A*c.^3))))-1/1728;
      coneq(260)=b'*((A*c.*(A*c.*(A*(A*c.^2)))))-1/5184;
      coneq(261)=b'*((A*c.*(A*c.*(A*(A*(A*c))))))-1/10368;
      coneq(262)=b'*((A*c.*(A*c.*(A*c.*(A*c)))))-1/3456;
      coneq(263)=b'*((A*c.*(A*c.^2.*(A*c.^2))))-1/1296;
      coneq(264)=b'*((A*c.*(A*c.^2.*(A*(A*c)))))-1/2592;
      coneq(265)=b'*((A*c.*(A*c.^3.*(A*c))))-1/864;
      coneq(266)=b'*((A*c.*(A*(A*c).*(A*c.^2))))-1/2592;
      coneq(267)=b'*((A*c.*(A*(A*c).*(A*(A*c)))))-1/5184;
      coneq(268)=b'*((A*c.*(A*c.*(A*c).*(A*c))))-1/1728;
      coneq(269)=b'*((A*c.^2.*(A*c.^4)))-1/360;
      coneq(270)=b'*((A*c.^2.*(A*(A*c.^3))))-1/1440;
      coneq(271)=b'*((A*c.^2.*(A*(A*(A*c.^2)))))-1/4320;
      coneq(272)=b'*((A*c.^2.*(A*(A*(A*(A*c))))))-1/8640;
      coneq(273)=b'*((A*c.^2.*(A*(A*c.*(A*c)))))-1/2880;
      coneq(274)=b'*((A*c.^2.*(A*c.*(A*c.^2))))-1/1080;
      coneq(275)=b'*((A*c.^2.*(A*c.*(A*(A*c)))))-1/2160;
      coneq(276)=b'*((A*c.^2.*(A*c.^2.*(A*c))))-1/720;
      coneq(277)=b'*((A*c.^2.*(A*(A*c).*(A*c))))-1/1440;
      coneq(278)=b'*((A*c.^3.*(A*c.^3)))-1/288;
      coneq(279)=b'*((A*c.^3.*(A*(A*c.^2))))-1/864;
      coneq(280)=b'*((A*c.^3.*(A*(A*(A*c)))))-1/1728;
      coneq(281)=b'*((A*c.^3.*(A*c.*(A*c))))-1/576;
      coneq(282)=b'*((A*c.^4.*(A*c.^2)))-1/216;
      coneq(283)=b'*((A*c.^4.*(A*(A*c))))-1/432;
      coneq(284)=b'*((A*c.^5.*(A*c)))-1/144;
      coneq(285)=b'*((A*(A*c.^2).*(A*c.^3)))-1/864;
      coneq(286)=b'*((A*(A*c.^2).*(A*(A*c.^2))))-1/2592;
      coneq(287)=b'*((A*(A*c.^2).*(A*(A*(A*c)))))-1/5184;
      coneq(288)=b'*((A*(A*c.^2).*(A*c.*(A*c))))-1/1728;
      coneq(289)=b'*((A*(A*(A*c)).*(A*c.^3)))-1/1728;
      coneq(290)=b'*((A*(A*(A*c)).*(A*(A*c.^2))))-1/5184;
      coneq(291)=b'*((A*(A*(A*c)).*(A*(A*(A*c)))))-1/10368;
      coneq(292)=b'*((A*(A*(A*c)).*(A*c.*(A*c))))-1/3456;
      coneq(293)=b'*((A*(A*c).*(A*c.^4)))-1/720;
      coneq(294)=b'*((A*(A*c).*(A*(A*c.^3))))-1/2880;
      coneq(295)=b'*((A*(A*c).*(A*(A*(A*c.^2)))))-1/8640;
      coneq(296)=b'*((A*(A*c).*(A*(A*(A*(A*c))))))-1/17280;
      coneq(297)=b'*((A*(A*c).*(A*(A*c.*(A*c)))))-1/5760;
      coneq(298)=b'*((A*(A*c).*(A*c.*(A*c.^2))))-1/2160;
      coneq(299)=b'*((A*(A*c).*(A*c.*(A*(A*c)))))-1/4320;
      coneq(300)=b'*((A*(A*c).*(A*c.^2.*(A*c))))-1/1440;
      coneq(301)=b'*((A*(A*c).*(A*(A*c).*(A*c))))-1/2880;
      coneq(302)=b'*((A*c.*(A*c.^2).*(A*c.^2)))-1/648;
      coneq(303)=b'*((A*c.*(A*c.^2).*(A*(A*c))))-1/1296;
      coneq(304)=b'*((A*c.*(A*(A*c)).*(A*(A*c))))-1/2592;
      coneq(305)=b'*((A*c.*(A*c).*(A*c.^3)))-1/576;
      coneq(306)=b'*((A*c.*(A*c).*(A*(A*c.^2))))-1/1728;
      coneq(307)=b'*((A*c.*(A*c).*(A*(A*(A*c)))))-1/3456;
      coneq(308)=b'*((A*c.*(A*c).*(A*c.*(A*c))))-1/1152;
      coneq(309)=b'*((A*c.^2.*(A*c).*(A*c.^2)))-1/432;
      coneq(310)=b'*((A*c.^2.*(A*c).*(A*(A*c))))-1/864;
      coneq(311)=b'*((A*c.^3.*(A*c).*(A*c)))-1/288;
      coneq(312)=b'*((A*(A*c).*(A*c).*(A*c.^2)))-1/864;
      coneq(313)=b'*((A*(A*c).*(A*c).*(A*(A*c))))-1/1728;
      coneq(314)=b'*((A*c.*(A*c).*(A*c).*(A*c)))-1/576;
      coneq(315)=b'*((A*c.^6).*c)-1/63;
      coneq(316)=b'*((A*(A*c.^5)).*c)-1/378;
      coneq(317)=b'*((A*(A*(A*c.^4))).*c)-1/1890;
      coneq(318)=b'*((A*(A*(A*(A*c.^3)))).*c)-1/7560;
      coneq(319)=b'*((A*(A*(A*(A*(A*c.^2))))).*c)-1/22680;
      coneq(320)=b'*((A*(A*(A*(A*(A*(A*c)))))).*c)-1/45360;
      coneq(321)=b'*((A*(A*(A*(A*c.*(A*c))))).*c)-1/15120;
      coneq(322)=b'*((A*(A*(A*c.*(A*c.^2)))).*c)-1/5670;
      coneq(323)=b'*((A*(A*(A*c.*(A*(A*c))))).*c)-1/11340;
      coneq(324)=b'*((A*(A*(A*c.^2.*(A*c)))).*c)-1/3780;
      coneq(325)=b'*((A*(A*(A*(A*c).*(A*c)))).*c)-1/7560;
      coneq(326)=b'*((A*(A*c.*(A*c.^3))).*c)-1/1512;
      coneq(327)=b'*((A*(A*c.*(A*(A*c.^2)))).*c)-1/4536;
      coneq(328)=b'*((A*(A*c.*(A*(A*(A*c))))).*c)-1/9072;
      coneq(329)=b'*((A*(A*c.*(A*c.*(A*c)))).*c)-1/3024;
      coneq(330)=b'*((A*(A*c.^2.*(A*c.^2))).*c)-1/1134;
      coneq(331)=b'*((A*(A*c.^2.*(A*(A*c)))).*c)-1/2268;
      coneq(332)=b'*((A*(A*c.^3.*(A*c))).*c)-1/756;
      coneq(333)=b'*((A*(A*(A*c).*(A*c.^2))).*c)-1/2268;
      coneq(334)=b'*((A*(A*(A*c).*(A*(A*c)))).*c)-1/4536;
      coneq(335)=b'*((A*(A*c.*(A*c).*(A*c))).*c)-1/1512;
      coneq(336)=b'*((A*c.*(A*c.^4)).*c)-1/315;
      coneq(337)=b'*((A*c.*(A*(A*c.^3))).*c)-1/1260;
      coneq(338)=b'*((A*c.*(A*(A*(A*c.^2)))).*c)-1/3780;
      coneq(339)=b'*((A*c.*(A*(A*(A*(A*c))))).*c)-1/7560;
      coneq(340)=b'*((A*c.*(A*(A*c.*(A*c)))).*c)-1/2520;
      coneq(341)=b'*((A*c.*(A*c.*(A*c.^2))).*c)-1/945;
      coneq(342)=b'*((A*c.*(A*c.*(A*(A*c)))).*c)-1/1890;
      coneq(343)=b'*((A*c.*(A*c.^2.*(A*c))).*c)-1/630;
      coneq(344)=b'*((A*c.*(A*(A*c).*(A*c))).*c)-1/1260;
      coneq(345)=b'*((A*c.^2.*(A*c.^3)).*c)-1/252;
      coneq(346)=b'*((A*c.^2.*(A*(A*c.^2))).*c)-1/756;
      coneq(347)=b'*((A*c.^2.*(A*(A*(A*c)))).*c)-1/1512;
      coneq(348)=b'*((A*c.^2.*(A*c.*(A*c))).*c)-1/504;
      coneq(349)=b'*((A*c.^3.*(A*c.^2)).*c)-1/189;
      coneq(350)=b'*((A*c.^3.*(A*(A*c))).*c)-1/378;
      coneq(351)=b'*((A*c.^4.*(A*c)).*c)-1/126;
      coneq(352)=b'*((A*(A*c.^2).*(A*c.^2)).*c)-1/567;
      coneq(353)=b'*((A*(A*c.^2).*(A*(A*c))).*c)-1/1134;
      coneq(354)=b'*((A*(A*(A*c)).*(A*(A*c))).*c)-1/2268;
      coneq(355)=b'*((A*(A*c).*(A*c.^3)).*c)-1/504;
      coneq(356)=b'*((A*(A*c).*(A*(A*c.^2))).*c)-1/1512;
      coneq(357)=b'*((A*(A*c).*(A*(A*(A*c)))).*c)-1/3024;
      coneq(358)=b'*((A*(A*c).*(A*c.*(A*c))).*c)-1/1008;
      coneq(359)=b'*((A*c.*(A*c).*(A*c.^2)).*c)-1/378;
      coneq(360)=b'*((A*c.*(A*c).*(A*(A*c))).*c)-1/756;
      coneq(361)=b'*((A*c.^2.*(A*c).*(A*c)).*c)-1/252;
      coneq(362)=b'*((A*(A*c).*(A*c).*(A*c)).*c)-1/504;
      coneq(363)=b'*((A*c.^5).*c.^2)-1/54;
      coneq(364)=b'*((A*(A*c.^4)).*c.^2)-1/270;
      coneq(365)=b'*((A*(A*(A*c.^3))).*c.^2)-1/1080;
      coneq(366)=b'*((A*(A*(A*(A*c.^2)))).*c.^2)-1/3240;
      coneq(367)=b'*((A*(A*(A*(A*(A*c))))).*c.^2)-1/6480;
      coneq(368)=b'*((A*(A*(A*c.*(A*c)))).*c.^2)-1/2160;
      coneq(369)=b'*((A*(A*c.*(A*c.^2))).*c.^2)-1/810;
      coneq(370)=b'*((A*(A*c.*(A*(A*c)))).*c.^2)-1/1620;
      coneq(371)=b'*((A*(A*c.^2.*(A*c))).*c.^2)-1/540;
      coneq(372)=b'*((A*(A*(A*c).*(A*c))).*c.^2)-1/1080;
      coneq(373)=b'*((A*c.*(A*c.^3)).*c.^2)-1/216;
      coneq(374)=b'*((A*c.*(A*(A*c.^2))).*c.^2)-1/648;
      coneq(375)=b'*((A*c.*(A*(A*(A*c)))).*c.^2)-1/1296;
      coneq(376)=b'*((A*c.*(A*c.*(A*c))).*c.^2)-1/432;
      coneq(377)=b'*((A*c.^2.*(A*c.^2)).*c.^2)-1/162;
      coneq(378)=b'*((A*c.^2.*(A*(A*c))).*c.^2)-1/324;
      coneq(379)=b'*((A*c.^3.*(A*c)).*c.^2)-1/108;
      coneq(380)=b'*((A*(A*c).*(A*c.^2)).*c.^2)-1/324;
      coneq(381)=b'*((A*(A*c).*(A*(A*c))).*c.^2)-1/648;
      coneq(382)=b'*((A*c.*(A*c).*(A*c)).*c.^2)-1/216;
      coneq(383)=b'*((A*c.^4).*c.^3)-1/45;
      coneq(384)=b'*((A*(A*c.^3)).*c.^3)-1/180;
      coneq(385)=b'*((A*(A*(A*c.^2))).*c.^3)-1/540;
      coneq(386)=b'*((A*(A*(A*(A*c)))).*c.^3)-1/1080;
      coneq(387)=b'*((A*(A*c.*(A*c))).*c.^3)-1/360;
      coneq(388)=b'*((A*c.*(A*c.^2)).*c.^3)-1/135;
      coneq(389)=b'*((A*c.*(A*(A*c))).*c.^3)-1/270;
      coneq(390)=b'*((A*c.^2.*(A*c)).*c.^3)-1/90;
      coneq(391)=b'*((A*(A*c).*(A*c)).*c.^3)-1/180;
      coneq(392)=b'*((A*c.^3).*c.^4)-1/36;
      coneq(393)=b'*((A*(A*c.^2)).*c.^4)-1/108;
      coneq(394)=b'*((A*(A*(A*c))).*c.^4)-1/216;
      coneq(395)=b'*((A*c.*(A*c)).*c.^4)-1/72;
      coneq(396)=b'*((A*c.^2).*c.^5)-1/27;
      coneq(397)=b'*((A*(A*c)).*c.^5)-1/54;
      coneq(398)=b'*((A*c).*c.^6)-1/18;
      coneq(399)=b'*((A*c.^3).*(A*c.^3))-1/144;
      coneq(400)=b'*((A*(A*c.^2)).*(A*c.^3))-1/432;
      coneq(401)=b'*((A*(A*(A*c))).*(A*c.^3))-1/864;
      coneq(402)=b'*((A*c.*(A*c)).*(A*c.^3))-1/288;
      coneq(403)=b'*((A*(A*c.^2)).*(A*(A*c.^2)))-1/1296;
      coneq(404)=b'*((A*(A*(A*c))).*(A*(A*c.^2)))-1/2592;
      coneq(405)=b'*((A*c.*(A*c)).*(A*(A*c.^2)))-1/864;
      coneq(406)=b'*((A*(A*(A*c))).*(A*(A*(A*c))))-1/5184;
      coneq(407)=b'*((A*c.*(A*c)).*(A*(A*(A*c))))-1/1728;
      coneq(408)=b'*((A*c.*(A*c)).*(A*c.*(A*c)))-1/576;
      coneq(409)=b'*((A*c.^4).*(A*c.^2))-1/135;
      coneq(410)=b'*((A*(A*c.^3)).*(A*c.^2))-1/540;
      coneq(411)=b'*((A*(A*(A*c.^2))).*(A*c.^2))-1/1620;
      coneq(412)=b'*((A*(A*(A*(A*c)))).*(A*c.^2))-1/3240;
      coneq(413)=b'*((A*(A*c.*(A*c))).*(A*c.^2))-1/1080;
      coneq(414)=b'*((A*c.*(A*c.^2)).*(A*c.^2))-1/405;
      coneq(415)=b'*((A*c.*(A*(A*c))).*(A*c.^2))-1/810;
      coneq(416)=b'*((A*c.^2.*(A*c)).*(A*c.^2))-1/270;
      coneq(417)=b'*((A*(A*c).*(A*c)).*(A*c.^2))-1/540;
      coneq(418)=b'*((A*c.^4).*(A*(A*c)))-1/270;
      coneq(419)=b'*((A*(A*c.^3)).*(A*(A*c)))-1/1080;
      coneq(420)=b'*((A*(A*(A*c.^2))).*(A*(A*c)))-1/3240;
      coneq(421)=b'*((A*(A*(A*(A*c)))).*(A*(A*c)))-1/6480;
      coneq(422)=b'*((A*(A*c.*(A*c))).*(A*(A*c)))-1/2160;
      coneq(423)=b'*((A*c.*(A*c.^2)).*(A*(A*c)))-1/810;
      coneq(424)=b'*((A*c.*(A*(A*c))).*(A*(A*c)))-1/1620;
      coneq(425)=b'*((A*c.^2.*(A*c)).*(A*(A*c)))-1/540;
      coneq(426)=b'*((A*(A*c).*(A*c)).*(A*(A*c)))-1/1080;
      coneq(427)=b'*((A*c.^5).*(A*c))-1/108;
      coneq(428)=b'*((A*(A*c.^4)).*(A*c))-1/540;
      coneq(429)=b'*((A*(A*(A*c.^3))).*(A*c))-1/2160;
      coneq(430)=b'*((A*(A*(A*(A*c.^2)))).*(A*c))-1/6480;
      coneq(431)=b'*((A*(A*(A*(A*(A*c))))).*(A*c))-1/12960;
      coneq(432)=b'*((A*(A*(A*c.*(A*c)))).*(A*c))-1/4320;
      coneq(433)=b'*((A*(A*c.*(A*c.^2))).*(A*c))-1/1620;
      coneq(434)=b'*((A*(A*c.*(A*(A*c)))).*(A*c))-1/3240;
      coneq(435)=b'*((A*(A*c.^2.*(A*c))).*(A*c))-1/1080;
      coneq(436)=b'*((A*(A*(A*c).*(A*c))).*(A*c))-1/2160;
      coneq(437)=b'*((A*c.*(A*c.^3)).*(A*c))-1/432;
      coneq(438)=b'*((A*c.*(A*(A*c.^2))).*(A*c))-1/1296;
      coneq(439)=b'*((A*c.*(A*(A*(A*c)))).*(A*c))-1/2592;
      coneq(440)=b'*((A*c.*(A*c.*(A*c))).*(A*c))-1/864;
      coneq(441)=b'*((A*c.^2.*(A*c.^2)).*(A*c))-1/324;
      coneq(442)=b'*((A*c.^2.*(A*(A*c))).*(A*c))-1/648;
      coneq(443)=b'*((A*c.^3.*(A*c)).*(A*c))-1/216;
      coneq(444)=b'*((A*(A*c).*(A*c.^2)).*(A*c))-1/648;
      coneq(445)=b'*((A*(A*c).*(A*(A*c))).*(A*c))-1/1296;
      coneq(446)=b'*((A*c.*(A*c).*(A*c)).*(A*c))-1/432;
      coneq(447)=b'*((A*c.^3).*(A*c.^2).*c)-1/108;
      coneq(448)=b'*((A*(A*c.^2)).*(A*c.^2).*c)-1/324;
      coneq(449)=b'*((A*(A*(A*c))).*(A*c.^2).*c)-1/648;
      coneq(450)=b'*((A*c.*(A*c)).*(A*c.^2).*c)-1/216;
      coneq(451)=b'*((A*c.^3).*(A*(A*c)).*c)-1/216;
      coneq(452)=b'*((A*(A*c.^2)).*(A*(A*c)).*c)-1/648;
      coneq(453)=b'*((A*(A*(A*c))).*(A*(A*c)).*c)-1/1296;
      coneq(454)=b'*((A*c.*(A*c)).*(A*(A*c)).*c)-1/432;
      coneq(455)=b'*((A*c.^4).*(A*c).*c)-1/90;
      coneq(456)=b'*((A*(A*c.^3)).*(A*c).*c)-1/360;
      coneq(457)=b'*((A*(A*(A*c.^2))).*(A*c).*c)-1/1080;
      coneq(458)=b'*((A*(A*(A*(A*c)))).*(A*c).*c)-1/2160;
      coneq(459)=b'*((A*(A*c.*(A*c))).*(A*c).*c)-1/720;
      coneq(460)=b'*((A*c.*(A*c.^2)).*(A*c).*c)-1/270;
      coneq(461)=b'*((A*c.*(A*(A*c))).*(A*c).*c)-1/540;
      coneq(462)=b'*((A*c.^2.*(A*c)).*(A*c).*c)-1/180;
      coneq(463)=b'*((A*(A*c).*(A*c)).*(A*c).*c)-1/360;
      coneq(464)=b'*((A*c.^2).*(A*c.^2).*c.^2)-1/81;
      coneq(465)=b'*((A*(A*c)).*(A*c.^2).*c.^2)-1/162;
      coneq(466)=b'*((A*(A*c)).*(A*(A*c)).*c.^2)-1/324;
      coneq(467)=b'*((A*c.^3).*(A*c).*c.^2)-1/72;
      coneq(468)=b'*((A*(A*c.^2)).*(A*c).*c.^2)-1/216;
      coneq(469)=b'*((A*(A*(A*c))).*(A*c).*c.^2)-1/432;
      coneq(470)=b'*((A*c.*(A*c)).*(A*c).*c.^2)-1/144;
      coneq(471)=b'*((A*c.^2).*(A*c).*c.^3)-1/54;
      coneq(472)=b'*((A*(A*c)).*(A*c).*c.^3)-1/108;
      coneq(473)=b'*((A*c).*(A*c).*c.^4)-1/36;
      coneq(474)=b'*((A*c.^3).*(A*c).*(A*c))-1/144;
      coneq(475)=b'*((A*(A*c.^2)).*(A*c).*(A*c))-1/432;
      coneq(476)=b'*((A*(A*(A*c))).*(A*c).*(A*c))-1/864;
      coneq(477)=b'*((A*c.*(A*c)).*(A*c).*(A*c))-1/288;
      coneq(478)=b'*((A*c.^2).*(A*c.^2).*(A*c))-1/162;
      coneq(479)=b'*((A*(A*c)).*(A*c.^2).*(A*c))-1/324;
      coneq(480)=b'*((A*(A*c)).*(A*(A*c)).*(A*c))-1/648;
      coneq(481)=b'*((A*c.^2).*(A*c).*(A*c).*c)-1/108;
      coneq(482)=b'*((A*(A*c)).*(A*c).*(A*c).*c)-1/216;
      coneq(483)=b'*((A*c).*(A*c).*(A*c).*c.^2)-1/72;
      coneq(484)=b'*((A*c).*(A*c).*(A*c).*(A*c))-1/144;
      coneq(485)=b'*(c.^8)-1/9;
    end
    if p>9
      disp('Order conditions for p>9 are not coded up yet');
    end

    coneq(end+1,1) = 0;
    coneq(2:end) = coneq(1:end-1);
    coneq(1) = sum(b) - 1;

    if ordernumFlag
        coneq = length(coneq);
    end
end