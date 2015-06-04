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

beta = 1;
   
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
        m = ceil(T/h);                % time points
        dt = T/m;                     % time steps
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
       
        % compute 1D matrices in hat basis
        e = ones(n, 1);
	    eshort = ones(n-2, 1);
	    Am = h/6 * spdiags([e [2; 4*eshort; 2] e], [-1 0 1], n, n); % mass matrix
	    As = 1/h * spdiags([-e [1; 2*eshort; 1] -e], [-1 0 1], n, n); % stiffness matrix
	    Ac = 1/2 * spdiags([-e [-1; 0*eshort; 1] e], [-1 0 1], n, n); % cross matrix
    
        % tensor products
        A = kron(Q(1,1)/2*As+mu(1)*Ac+r*Am,Am) + kron((-Q(1,2)*Ac+mu(2)*Am),Ac) + Q(2,2)/2*kron(Am,As);
        Am = kron(Am,Am); 
                
        % projection for payoff
        fct = @(x1,x2) (max(0,a*exp(x1)-b*exp(x2)));
        f = rhs2d(x,fct);
        u0 = Am\f;                        
                                                      
        %------------------------------------------------------------------
        %  Solver
        %------------------------------------------------------------------

        theta = 0.5;
            
        % graded time mesh
        t = (0:dt:T).^beta;

        % full grid
        u = u0; maxiter = 0;
        for i=1:m 
            B1 = Am-(1-theta)*(t(i+1)-t(i))*A;
            B2 = Am+theta*(t(i+1)-t(i))*A;
            [u,flag,res,iter] = gmres(B2,B1*u,[],1.0e-7,min([200 size(B2,1)]),[],[],u);
            maxiter = max([maxiter iter]);
            %u=B2/B1*u;
        end

        %fprintf('Maximum of %2d GMRES iterations\n',maxiter)
        u = reshape(u,n,n);

        % area of interest
        K = 10;
        S1 = reshape(S1,n,n); 
        S2 = reshape(S2,n,n);
        u0 = max(a*S1-b*S2,0);
                
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
        u = u(I2,I1); u0 = max(a*S1-b*S2,0);
        % compute L_infty-error 
        exact = bs_exchange([S1(:) S2(:)],T,sigma,rho,a,b); 
        errinf(j) = max(exact - u(:));

    end

    % elapsed time
    toc 

%--------------------------------------------------------------------------
%  Postprocessing
%--------------------------------------------------------------------------

% plot error with respect to h
H = 2.^(-L-1)';
figure(1)
h = axes;
loglog(H,err,'bx-');
hold on
loglog(H,errinf,'rx-');
hold on
p = polyfit(log(H),log(err),1);
loc = ginput(1);
mySlope(h,loc,0,p(1),'k+-')
legend('Exchange option','Location','NorthWest')
title('Error with respect to h')

set(h,'FontSize',14);
axis on
xlabel('h')
ylabel('Linf-Error')
hold off

% plot error with respect to N
% compare with (8.24) in the book
N = (2.^(L'+1)).^2;
figure(2); clf;
h = axes;
loglog(N,err,'bx-');
hold on
loglog(N,errinf,'rx-');
hold on
p = polyfit(log(N),log(err),1);
loc = ginput(1);
mySlope(h,loc,1,p(1),'k+-')
legend('Exchange option','Location','NorthEast')
title('Error with respect to N')

set(h,'FontSize',14);
axis on
xlabel('N')
ylabel('Linf-Error')
hold off

exact = bs_exchange([S1(:) S2(:)],T,sigma,rho,a,b);
exact = reshape(exact,length(S1),length(S2));

% plot option price
figure(3)
subplot(2,1,1)
surf(S1,S2,u)
hold on
mesh(S1,S2,u0)
hold off
box on
grid on

set(h,'FontSize',14)
axis on
xlabel('s_1')
ylabel('s_2')
zlabel('FE Option price','FontSize',14)

subplot(2,1,2)
surf(S1,S2,exact)
hold on
mesh(S1,S2,u0)
hold off
box on
grid on

set(h,'FontSize',14)
axis on
xlabel('s_1')
ylabel('s_2')
zlabel('Exact Option price','FontSize',14)
end






