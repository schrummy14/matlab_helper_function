% RK23 max Eff
function [res,res1,res2,A,b,bp,c] = rk23Finder(x,normVal,pow)
 % Possible x value
%  nx =[
%     8469940   151404220    35793988
%     ];
% 
% dx = [
%     13215339   915854223   347153397
%     ];

 
bp3 = x(3);
b3  = x(2);
bp2 = x(1);

c = [
    0,...
    (2*(-1+3*b3))/(3*(-1+2*b3)),...
    1,...
    1];

b = [
    (-1+4*b3)/(4*(-1+3*b3)),...
    -(3*(1-4*b3+4*b3^2)/(4*(-1+3*b3))),...
    b3,...
    0];

bp =[
    (-3+6*b3+2*bp2)/(6*(-1+2*b3)),...
    bp2,...
    bp3,...
    (-3+6*b3+4*bp2-12*b3*bp2+6*bp3-12*b3*bp3)/(6*(-1+2*b3))];

a32 = ((-1+2*b3)/(4*b3*(-1+3*b3)));

A = [
    0         0    0    0
    c(2)      0    0    0
    c(3)-a32  a32  0    0
    b(1)      b(2) b(3) b(4) 
    ];  

res1 = oc_butcher(A,b,c,pow);
res2 = oc_butcher(A,bp,c,pow);

res = norm([res1;res2],normVal);