% Input: u0: payoff function for different values of S1 and S2
function u = PDESolver(x, n, T, h, Q, mu, r, f,nb)

beta = 1;

        m = ceil(T/h);                % time points
        dt = T/m;                     % time steps
        
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
        %f = rhs2d(x,u0);
        if (nb == 1)
        u0 = Am\f;   
        else
        u0 = f;
        end
                                                      
        %------------------------------------------------------------------
        %  Solver
        %------------------------------------------------------------------

        theta = 0.5;
            
        % graded time mesh
        t = (0:dt:T).^beta;

        % full grid
        u = u0; 
        maxiter = 0;
        for i=1:m
            B1 = Am-(1-theta)*(t(i+1)-t(i))*A;
            B2 = Am+theta*(t(i+1)-t(i))*A;          
            [u,flag,res,iter] = gmres(B2,B1*u,[],1.0e-7,min([200 size(B2,1)]),[],[],u);
            maxiter = max([maxiter iter]);
            %u=B2/B1*u;
        end

        %fprintf('Maximum of %2d GMRES iterations\n',maxiter)
       % u = reshape(u,n,n);

u;
end






