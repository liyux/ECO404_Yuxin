clear all
close all

%keep constants as scalars, make into vectors only when needed

%coefficients
%coef matrix is 1x6
coef.constant=0.5184;
coef.lawnsize=0.0231;
coef.weather=1.5945;
coef.bath=0.1336;
coef.hsize=0.0046;
coef.billing=0.034;
coef.p1=-1.8989;
coef.yd=0.1782;
coefmatrix=[coef.constant coef.lawnsize coef.weather coef.bath coef.hsize coef.billing];

%fixed variables
numhh=10;
constant=1;
billing=30;
weather=0.464;
p11=1.2;
p12=1.4;
x11=12;
d1=0;
d2=0.001*(p12-p11)*x11;

correlationmatrix=[1 0.6 0.7 0.7 ; 0.6 1 0.5 0.5 ; 0.7 0.5 1 0.9 ; 0.7 0.5 0.9 1];

%household characteristics
%to keep values inside a range, generate more than needed, remove those
%outside, and keep top 1000 numbers. csubset=condd1*find(condd1>x11) is subset,
%keep csubset(1:1000,1)

income=normrnd(1.156,0.485,[1.5*numhh,1]); %$1000
isubset=income(find((income>=0.137)&(income<=2.967)));
income = isubset(1:numhh);
hsize=normrnd(18.337,5.288,[1.5*numhh,1]); %100sqft
hsubset=hsize(find((hsize>=4.44)&(hsize<=36.11)));
hsize = hsubset(1:numhh);
lawnsize=normrnd(9.877,3.383,[1.5*numhh,1]); %1000sqft
lsubset=lawnsize(find((lawnsize>=4.61)&(lawnsize<=25.96)));
lawnsize = lsubset(1:numhh);
bath=normrnd(1.635,0.52,[1.5*numhh,1]);
bsubset=bath(find((bath>=1)&(bath<=3)));
bath = bsubset(1:numhh);

correlatedhh=[income lawnsize bath hsize]*chol(correlationmatrix);
correlatedhh=sortrows(correlatedhh,1);
income=correlatedhh(1:numhh,1);
hhmatrix=[constant*ones(numhh,1) correlatedhh(1:numhh,2) weather*ones(numhh,1) correlatedhh(1:numhh,3) correlatedhh(1:numhh,4) billing*ones(numhh,1)];
%hhmatrix is hhsizex6


%calculations to get IBP demand
d.hh=exp(hhmatrix*coefmatrix'); %size numhhx1
d.p11=p11.^coef.p1;
d.p12=p12.^coef.p1;
d.y1=(income+d1).^coef.yd;
d.y2=(income+d2).^coef.yd;

condd1=d.hh.*d.p11.*d.y1; %size numhhx1
condd2=d.hh.*d.p12.*d.y2; %size numhhx1
condd = [condd1 ones(numhh,1)*x11 condd2];

% changing p1
p11exp = 1.1:0.01:1.3;
numIterp1 = numel(p11exp);
qtydIBPp1 = zeros(numhh,numIterp1);
expIBPp1= zeros(numhh,numIterp1);

for i = 1:numIterp1;
    
    p11i = p11exp(i);
    %redefining all variables containing p11
    d2 = 0.001*(p12-p11exp(i))*x11;
    pricematrix = [p11exp(i)*ones(numhh,1) p11exp(i)*ones(numhh,1) p12*ones(numhh,1)];
    d.p11 = p11exp(i).^coef.p1;
    d.y2 = (income+d2).^coef.yd;
    condd1 = d.hh.*d.p11.*d.y1;
    condd2 = d.hh.*d.p12.*d.y2;
    condd = [condd1 ones(numhh,1)*x11 condd2];
    
    % creating a matrix of binary variables that has the value 1 if the
    %household belongs in that block, and 0 otherwise. 
    
    checkCondd1 = condd1<x11;
    checkCondd2 = condd2>x11;
    checkAtKink = ones(numhh,1) - checkCondd1 - checkCondd2;
    hhBlock = [checkCondd1 checkAtKink checkCondd2];
    qtydIBPp1(:,i) = sum(hhBlock.*condd,2);
    
    % finding total expenditure per household by multiplying prices with
    % quantity demanded in each block. Total revenue is obtained by summing
    % the column of household expenditures.
    
    expIBPp1(:,i) = sum(hhBlock.*condd.*pricematrix,2);
    
end

%changing x11
x11exp = 10:1:20;
numIterx1 = numel(x11exp);
qtydIBPx1 = zeros(numhh,numIterx1);
expIBPx1= zeros(numhh,numIterx1);

%for p11 = 1.1:0.01:1.3
for j = 1:numIterx1;
    
    x11i = x11exp(j);
    %redefining all variables containing x11
    d2 = 0.001*(p12-p11)*x11exp(j);
    d.y2=(income+d2).^coef.yd;
    condd2=d.hh.*d.p12.*d.y2;
    condd = [condd1 ones(numhh,1)*x11exp(j) condd2];
    
    % creating a matrix of binary variables that has the value 1 if the
    %household belongs in that block, and 0 otherwise. 
    
    checkCondd1 = condd1<x11exp(j);
    checkCondd2 = condd2>x11exp(j);
    checkAtKink = ones(numhh,1) - checkCondd1 - checkCondd2;
    hhBlock = [checkCondd1 checkAtKink checkCondd2];
    qtydIBPx1(:,i) = sum(hhBlock.*condd,2);
    
    % finding total expenditure per household by multiplying prices with
    % quantity demanded in each block. Total revenue is obtained by summing
    % the column of household expenditures.
    
    expIBPx1(:,i) = sum(hhBlock.*condd.*pricematrix,2);
    
end

%move 2 things: to make same revenue at this price, where does the break
%point have to be?

%if exist ('hhmatrix.mat', 'file')
%load hhmatrix
%else (that chunk about generating, add in code that saves the resulting matrix (save ('hhmatrix','matrix'))

%quantile(vector, #breaks)
%redefine vector to only include those in a range
%sort vector and divide it

%what does the 0.5 /1 mean
%look for papers saying low income like IBP and why.