% MAIN_BSEXCHANGE computes convergence rate for the exchange option in 
% the BS market model 

function main
clear all;
close all;
tic

%--------------------------------------------------------------------------
%  Set Parameters 
%--------------------------------------------------------------------------
L = 5;                       % number of levels
R = 2;                       % domain (-R,R)^2
TC = [4,3,2,1];              % Callable dates (yearly)
%TCoupon = [0.5,1,1.5,2,2.5]; % Coupon dates (every 6 months)
T = 5;                       % maturity
a = 1;                       % constant in payoff
b = 2;                       % constant in payoff
r = 0.00;                   % interest rates
sigma = [0.4;0.2];           % volatilities
rho = 0.2;                   % correlation
K = 1;                       % strike
N = 1;                       % Notional
B = 0.7;                     % barrier level
C = 0.0;                     % Fixed coupon

Q = zeros(2,2);
Q(1,1) = sigma(1)^2;
Q(1,2) = sigma(1)*sigma(2)*rho;
Q(2,1) = Q(1,2);
Q(2,2) = sigma(2)^2;
mu = [Q(1,1)/2; Q(2,2)] -r;

fprintf('Level L = %d\n',L)
%------------------------------------------------------------------
%  Discretization
%------------------------------------------------------------------
                
% 1d mesh using barrier
nb = 2.^(L+1)-1;                 % number of nodes using barrier
h = (R - log(B))/(nb+1);         % mesh size
xb = linspace(log(B),R,nb+2)';   % mesh nodes
dofb = 2:nb+1;                   % degree of freedoms
        
% with barrier
eb = ones(nb,1);
XB1 = (xb(dofb)*eb')';
XB2 = xb(dofb)*eb';
SB1 = exp(XB1);
SB2 = exp(XB2);
        
IB1 = (abs(SB1(1,:)-K) < K/2);
IB2 = (abs(SB2(:,2)-K) < K/2);
SB1 = SB1(IB2,IB1);
SB2 = SB2(IB2,IB1);

% 1d mesh not using barrier
n = 2.^(L+2)-1;
Start = R-(n+1)*h;
x = linspace(Start,R,n+2)';
dof = 2:n+1;               
        
% 2d mesh
% axiparallell mesh
e = ones(n,1);
X1 = (x(dof)*e')';
X2 = x(dof)*e';
S1 = exp(X1);
S2 = exp(X2);

% area of interest
S1 = reshape(S1,n,n); 
S2 = reshape(S2,n,n);

I1 = (abs(S1(1,:)-K) < K/2);
I2 = (abs(S2(:,2)-K) < K/2);
S1 = S1(I1,I1);
S2 = S2(I1,I1);

%------------------------------------------------------------------
%  RHS
%------------------------------------------------------------------

fct = @(x1,x2) (-max(K - min(exp(x1),exp(x2)),0));
f = rhs2d(x,fct);
fb = rhs2d(xb,fct);
        
% Optional:
% TODO: Extend to 3 dimentions!
% TODO: Extend to stochastic volatility!
% TODO: Extend to non-constant correlation!
% TODO: Extend to include correlation between bond and option part!
% TODO: Extend to include stochastic interest rates!!!
% TODO: Extend to include counterparty credit risk!!!!

% First callable date

% put option: in case barrier has already been hit


% u1 is the do-put option
% we are looking for the di-put option

% If we are at the last callable date and keep the BRC until
% maturity we get bond - down-in put (if barrier has not been hit)
% or bond - put (if barrier has already been hit)

% If we call we get coupon and notional now


% => Value of BRC at last callable date!!

% If we are at the 2nd last callable date and barrier has been hit:
% Keep until next callable date:
% 
% Get value of callable at next callable date!
% U


% Barrier has not been hit:
%
% if called: notional + coupons now
%

InitialCondHit = f; % Payoff of option in case the barrier has been hit at maturity
InitialCondHitB = fb;
InitialCondNoHit = fb*0;        
ind = (n-nb+1):n;
dt = T - TC(1);
disc = exp(-r*dt);

% Further callable dates
firstTime = 1;
nbOfCallableDates = length(TC);
for (i = 1:nbOfCallableDates)        

    % Vanilla put (Hit)
    p = PDESolver(x, n, dt, h, Q, mu, r, InitialCondHit, firstTime);
    p = reshape(p,n,n);

    % Down-and-out put (Hit)
    do = PDESolver(xb, nb, dt, h, Q, mu, r, InitialCondHitB, firstTime);
    do = reshape(do,nb,nb);
    di = p;
    
    % Down-and-in put (Hit)
    di(ind,ind) = di(ind,ind) - do;

    % Down-and-out put (No Hit)
    do = zeros(n,n);            
    don = PDESolver(xb, nb, dt, h, Q, mu, r, InitialCondNoHit, firstTime);
    do(ind,ind) = reshape(don,nb,nb);
 
    % Not yet hit:
    % In case we hit barrier during this time period: already hit for next
    % time period
    % In case we don't hit barrier during this time period: not yet hit for
    % next time period
    % Down-and-in  (Hit) + Down-and-out (No hit)
    BRCnohitnocall = C + disc*(N+C) + di + do;
    BRCnohit = min(BRCnohitnocall, C+N);

    % Vanilla (Hit)
    BRChitnocall = C + p + disc*(N+C);
    BRChit = min(BRChitnocall, C + N);

    Diffhit = (BRChit - BRChitnocall);
    Diffnohit = BRCnohit - BRCnohitnocall; % More often different from 0            

    plotStuff(i,S1, S2, I1,I2,BRCnohitnocall, BRCnohit, BRChitnocall, BRChit, Diffnohit, Diffhit) 

    firstTime = 2;            
    % What you get at next time step in case barrier has already been hit
    InitialCondHit = BRChit - (N + C);
    InitialCondHitB = BRChit(ind,ind) - (N + C);
    % What BRC is worth at next time step in case barrier has not yet been
    % hit
    InitialCondNoHit = BRCnohit(ind,ind) - (N + C);
    if (i < nbOfCallableDates)
      dt = TC(i) - TC(i+1);
      disc = exp(-r*dt);
    end  
end

% elapsed time
toc 

%--------------------------------------------------------------------------
%  Postprocessing
%--------------------------------------------------------------------------



end
