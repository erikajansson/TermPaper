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
TC = [1,2];                  % Callable dates (yearly)
TCoupon = [0.5,1,1.5,2,2.5]; % Coupon dates (every 6 months)
T = 3;                       % maturity
a = 1;                       % constant in payoff
b = 2;                       % constant in payoff
r = 0;                       % interest rates
sigma = [0.4;0.1];           % volatilities
rho = -0.6;                  % correlation
K = 1;                       % strike
Notional = 1;
B = 0.8;                     % barrier level
Coupon = 0.05;               % Fixed coupon

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
        
        % with barrier
        eb = ones(nb,1);
        XB1 = (xb(dofb)*eb')';
        XB2 = xb(dofb)*eb';
        SB1 = exp(XB1);
	    SB2 = exp(XB2);
        
        %------------------------------------------------------------------
        %  Compute Stiffness Matrix
        %------------------------------------------------------------------

        fct = @(x1,x2) (max(K - min(exp(x1),exp(x2)),0));
        f = rhs2d(x,fct);
        fb = rhs2d(xb,fct);
        
        % Optional:
        % TODO: Extend to 3 dimentions!
        % TODO: Extend to stochastic volatility!
        % TODO: Extend to non-constant correlation!
        % TODO: Extend to include correlation between bond and option part!
        % TODO: Extend to include stochastic interest rates!!!
        % TODO: Extend to include counterparty credit risk!!!!   
        
        % TODO: Input: call-dates,
        
        % TODO: Add barrier !!!    
        % TODO: What should the barrier value be?        
              
        % First callable date
        
        pDO = PDESolver(xb, nb, TC(1), h, Q, mu, r, fb,1);
        pDO = reshape(pDO,nb,nb);
        p = PDESolver(x, n, TC(1), h, Q, mu, r, f,1);
        p = reshape(p,n,n);
        pDI = p;
        ind = (n-nb+1):n;
        size(pDI(ind,ind))
        size(pDO)
        pDI(ind,ind) = pDI(ind,ind) - pDO;
%       p - pDI
                
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
        
        % Further callable dates
%        u2 = u1;
%        nbOfCallableDates = length(TC);
%        for (i = 2:nbOfCallableDates)        
%          u2 = PDESolver(x, n, TC(i-1) - TC(i), h, Q, mu, r, u2, 2);           
%        end
        
        % TODO: Add bond part !!!!
        
        % TODO: Time between dates might be different !!!
        
%        u = PDESolver(x, n, T - TC(nbOfCallableDates), h, Q, mu, r, u2,2);        

%        u = reshape(u,n,n);
%        u1 = reshape(u1,n,n);
        % area of interest
        S1 = reshape(S1,n,n); 
        S2 = reshape(S2,n,n);
        u0 = max(K - min(S1,S2),0);
                
        s = 1;
                
        s1 = find(S1(1,:) == 1);
        s2 = find(S2(:,2) == 1);
        
        % exact = bs_exchange([s s],T,sigma,rho,a,b);
        % compute error at (1,1)
        % err(j) = exact - u(s1,s2);
                
        I1 = (abs(S1(1,:)-K) < K/2);
        I2 = (abs(S2(:,2)-K) < K/2);
        S1 = S1(I2,I1); 
        S2 = S2(I2,I1); 
%        u = u(I2,I1); 
        IB1 = (abs(SB1(1,:)-K) < K/2);
        IB2 = (abs(SB2(:,2)-K) < K/2);
        SB1 = SB1(IB2,IB1);
        SB2 = SB2(IB2,IB1);

    % elapsed time
    toc 

%--------------------------------------------------------------------------
%  Postprocessing
%--------------------------------------------------------------------------

% plot option price
figure(3)
%subplot(2,1,1)
%pDO = reshape(pDO,nb,nb);
surf(S1,S2,pDI(I1,I2))
hold on
%p = reshape(p,n,n);
%figure(4)

mesh(S1,S2,p(I1,I2))
hold off
box on
grid on

%set(h,'FontSize',14)
axis on
xlabel('s_1')
ylabel('s_2')
zlabel('FE Option price','FontSize',14)

end
