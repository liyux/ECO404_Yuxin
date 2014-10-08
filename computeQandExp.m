%calculate quantity demanded and expenditure within a loop

    %redefining all variables containing p11
    d2 = 0.001*(p12i-p11i)*x11i;
    pricematrix = [p11i*ones(numhh,1) p11i*ones(numhh,1) p12i*ones(numhh,1)];
    d.p11i = p11i.^coef.p1;
    d.p12i=p12i.^coef.p1;
    d.y2=(income+d2).^coef.yd;
    condd1 = d.hh.*d.p11i.*d.y1;
    condd2 = d.hh.*d.p12i.*d.y2;
    condd = [condd1 ones(numhh,1)*x11i condd2];
    
    % creating a matrix of binary variables that has the value 1 if the
    %household belongs in that block, and 0 otherwise. 
    
    checkCondd1 = condd1<x11i;
    checkCondd2 = condd2>x11i;
    checkAtKink = ones(numhh,1) - checkCondd1 - checkCondd2;
    hhBlock = [checkCondd1 checkAtKink checkCondd2];
    qtydIBPi(:,i) = sum(hhBlock.*condd,2);
    
    % finding total expenditure per household by multiplying prices with
    % quantity demanded in each block. Total revenue is obtained by summing
    % the column of household expenditures.
    
    expIBPi(:,i) = sum(hhBlock.*condd.*pricematrix,2);
