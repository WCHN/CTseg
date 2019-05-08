function scl = scl_from_bf(chan)
C   = numel(chan);
scl = zeros(1,C);
for c=1:C
    b1 = chan(c).B1(1,1);
    b2 = chan(c).B2(1,1);
    b3 = chan(c).B3(1,1);
    t1 = chan(c).T(1,1,1);
    
    scl(c) = exp(b1*b2*b3*t1);
end
%==========================================================================