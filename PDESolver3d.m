% Input: 
% u0: payoff function for different values of S1,S2 and S3
% h: grid size in space
% nb == 1: use projection for payoff!!!
function u = PDESolver3d(x, n, T, h, Q, mu, r, f,nb)

beta = 1;

m = ceil(T/h);                % time points
dt = T/m;                     % time steps
dof = 2:n-1;
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
MM = kron(Am,Am);
BM = kron(Ac,Am);
MB = kron(Am,Ac);
M1 = mu(3)*MM - Q(1,3)*BM - Q(2,3)*MB;
M2 = r*MM+mu(2)*MB+mu(1)*BM + 1/2*Q(1,1)*kron(As,Am)+1/2*Q(2,2)*kron(Am,As)-Q(1,2)*kron(Ac,Ac);

A = 1/2*Q(3,3)*kron(MM,As) + kron(M1,Ac) + kron(M2, Am);
Am = kron(MM, Am);

% projection for payoff
%f = rhs2d(x,u0);
size(Am)
f = f(:);
size(f)
if (nb == 1)
%  u0 = Am\f;   
  [u0,flag,res,iter] = gmres(Am,f,[],1.0e-7,min([200 size(Am,1)]),[],[],[]);
else
  u0 = f;
end

fprintf('Finished computing initial cond.\n')

%------------------------------------------------------------------
%  Solver
%------------------------------------------------------------------

theta = 0.5;

% graded time mesh
t = (0:dt:T).^beta;

bc = ones(n,n,n);
bc(1,1,:) = 0;
bc(n,n,:) = 0;
bc(1,n,:) = 0;
bc(n,1,:) = 0;            
bc(:,1,1) = 0;
bc(:,n,n) = 0;             
bc(:,1,n) = 0;
bc(:,n,1) = 0;
BCI = bc(:) == 0;

% full grid
u = u0; 
maxiter = 0;
tic
for i=1:m            
    B1 = Am-(1-theta)*(t(i+1)-t(i))*A;
    B2 = Am+theta*(t(i+1)-t(i))*A;   
    u(BCI) = 0;
    [u,flag,res,iter] = gmres(B2,B1*u,[],1.0e-7,min([200 size(B2,1)]),[],[],u);
    maxiter = max([maxiter iter]);
    %u=B2/B1*u;
end
toc
fprintf('Maximum of %2d GMRES iterations\n',maxiter)
% u = reshape(u,n,n);

u;
end






