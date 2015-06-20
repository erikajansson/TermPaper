% MAIN_BSEXCHANGE computes convergence rate for the exchange option in 
% the BS market model 

function main
clear all;
close all;
tic

%--------------------------------------------------------------------------
%  Set Parameters 
%--------------------------------------------------------------------------
L = 4;                        % number of levels
R = 2;                        % domain (-R,R)^2
TC = [4,3,2,1];               % Callable dates (yearly)
%TCoupon = [0.5,1,1.5,2,2.5]; % Coupon dates (every 6 months)
T = 5;                        % maturity
a = 1;                        % constant in payoff
b = 2;                        % constant in payoff
r = 0.00;                     % interest rates
sigma = [0.4;0.2;0.3];        % volatilities
rho = [0.2,0.8;0.2,0.5];      % correlation
K = 1;                        % strike
N = 1;                        % Notional
B = 0.7;                      % barrier level
C = 0.0;                      % Fixed coupon

Q = zeros(3,3);
Q(1,1) = sigma(1)^2;
Q(1,2) = sigma(1)*sigma(2)*rho(1,1);
Q(1,3) = sigma(1)*sigma(3)*rho(1,2);
Q(2,1) = Q(1,2);
Q(2,2) = sigma(2)^2;
Q(2,3) = sigma(2)*sigma(3)*rho(2,2);
Q(3,1) = Q(1,3);
Q(3,2) = Q(2,3);
Q(3,3) = sigma(3)^2;
mu = [Q(1,1); Q(2,2);Q(3,3)]/2 -r;

fprintf('Level L = %d\n',L)
%------------------------------------------------------------------
%  Discretization
%------------------------------------------------------------------
                
% 1d mesh using barrier
nb = 2.^(L+1)-1;                 % number of nodes using barrier
h = (R - log(B))/(nb+1);         % mesh size
xb = linspace(log(B),R,nb+2)';   % mesh nodes
dofb = 2:nb+1;                   % degree of freedoms
        
% 1d mesh not using barrier
n = 2.^(L+2)-1;
Start = R-(n+1)*h;
x = linspace(Start,R,n+2)';
dof = 2:n+1;       

I = find(x >= 0,1,'first');

%------------------------------------------------------------------
%  RHS
%-----------------------------------------------------------------

fct = @(x1,x2,x3) (-max(K - min(min(exp(x1),exp(x2)),exp(x3)),0));
fprintf('Evaluate f')
f = rhs3d(x,fct);
fprintf('Evaluate f for barrier')
f = reshape(f,n,n,n);
ind = (n-nb+1):n;
fb = f(ind,ind,ind);
%rhs3d(xb,fct);
fprintf('F Evaluated')
        
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
dt = T - TC(1);
disc = exp(-r*dt);

% Further callable dates
firstTime = 1;
nbOfCallableDates = length(TC);
for (i = 1:nbOfCallableDates)
    
    fprintf('For date nb %d\n', i)

    % Vanilla put (Hit)
    p = PDESolver3d(x, n, dt, h, Q, mu, r, InitialCondHit, firstTime);
    p = reshape(p,n,n,n);

    % Down-and-out put (Hit)
    do = PDESolver3d(xb, nb, dt, h, Q, mu, r, InitialCondHitB, firstTime);
    do = reshape(do,nb,nb,nb);
    di = p;
    
    % Down-and-in put (Hit)
    di(ind,ind,ind) = di(ind,ind,ind) - do;

    % Down-and-out put (No Hit)
    do = zeros(n,n);            
    don = PDESolver3d(xb, nb, dt, h, Q, mu, r, InitialCondNoHit, firstTime);
    do(ind,ind,ind) = reshape(don,nb,nb,nb);
 
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
   
    fprintf('Barrier already hit time step %d: %f\n',i,BRChit(I,I,I))
    fprintf('Barrier not yet hit time step %d: %f\n',i,BRCnohit(I,I,I))    
    
%    plotStuff(i,S1, S2, I1,I2,BRCnohitnocall, BRCnohit, BRChitnocall, BRChit, Diffnohit, Diffhit) 

    firstTime = 2;            
    % What you get at next time step in case barrier has already been hit
    InitialCondHit = BRChit - (N + C);
    InitialCondHitB = BRChit(ind,ind,ind) - (N + C);
    % What BRC is worth at next time step in case barrier has not yet been
    % hit
    InitialCondNoHit = BRCnohit(ind,ind,ind) - (N + C);
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
