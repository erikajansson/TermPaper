% MAIN_BSEXCHANGE computes convergence rate for the exchange option in 
% the BS market model 

function main_bsexchange
clear all;
close all;
tic

%--------------------------------------------------------------------------
%  Set Parameters 
%--------------------------------------------------------------------------

L = 4:6;                    % number of levels  
R = 5;                      % domain (-R,R)^2
T = 1;                      % maturity
a = 1;                      % constant in payoff
b = 2;                      % constant in payoff
r = 0;                      % interest rate 
sigma = [0.4;0.1];          % volatilities
rho = -0.6;                 % correlation
Q = zeros(2,2);
Q(1,1) = sigma(1)^2;
Q(1,2) = sigma(1)*sigma(2)*rho;
Q(2,1) = Q(1,2);
Q(2,2) = sigma(2)^2;
mu = [Q(1,1)/2; Q(2,2)] -r;
K = 1;
   
    % loop over levels
    err = zeros(length(L),1);
    errinf = zeros(length(L),1);
    for j=1:length(L)
    
        fprintf('Level L = %d\n',L(j))
    
        %------------------------------------------------------------------
        %  Discretization
        %------------------------------------------------------------------

        % 1d mesh
        n = 2.^(L(j)+1)-1;            % number of nodes
        h = 2*R/(n+1);                % mesh size
        x = linspace(-R,R,n+2)';      % mesh nodes        
        dof = 2:n+1;                  % degree of freedoms   

        % 2d mesh
        % axiparallell mesh
        e = ones(n,1);
        X1 = (x(dof)*e')';    
        X2 = x(dof)*e';
        S1 = exp(X1);
	    S2 = exp(X2);

        %------------------------------------------------------------------
        %  Compute Stiffness Matrix
        %------------------------------------------------------------------

        fct = @(x1,x2) (max(K - min(exp(x1),exp(x2)),0));
        %fct = @(x1,x2) (max(0,a*exp(x1)-b*exp(x2)));
        %u0 = fct(X1,X2);
        f = rhs2d(x,fct);          
        
        % compute 1D matrices in hat basis
        u = PDESolver(x, n, T, h, Q, mu, r, f);       

        % area of interest
        K = 10;
        S1 = reshape(S1,n,n); 
        S2 = reshape(S2,n,n);
        u0 = max(K - min(S1,S2),0);
        %fct(x,x);
        %size(u0)
%        u0 = reshape(u0,n,n);
%        max(a*S1-b*S2,0);
                
        s = 1;
                
        s1 = find(S1(1,:) == 1);
        s2 = find(S2(:,2) == 1);
        
        exact = bs_exchange([s s],T,sigma,rho,a,b);
        % compute error at (1,1)
        err(j) = exact - u(s1,s2);
                
        I1 = (abs(S1(1,:)-K) < K/2);
        I2 = (abs(S2(:,2)-K) < K/2);
        S1 = S1(I2,I1); 
        S2 = S2(I2,I1); 
        u = u(I2,I1); 
        %u0 = max(a*S1-b*S2,0);
        % compute L_infty-error 
        exact = bs_exchange([S1(:) S2(:)],T,sigma,rho,a,b); 
        errinf(j) = max(exact - u(:));

    end

    % elapsed time
    toc 

%--------------------------------------------------------------------------
%  Postprocessing
%--------------------------------------------------------------------------


% plot option price
figure(3)
%subplot(2,1,1)
surf(S1,S2,u)
hold on
mesh(S1,S2,u0(I1,I2))
hold off
box on
grid on

%set(h,'FontSize',14)
axis on
xlabel('s_1')
ylabel('s_2')
zlabel('FE Option price','FontSize',14)

end






