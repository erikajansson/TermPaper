%  RHS2D computes the right hand side / load vector (f,v) in
%     two dimension using 2D hat basis functions
%
%     L = rhs2d(vertices,FHandle,varargin)
%   
%     Input:  vertices  ... one-dimensional vertices
%             FHandle   ... function handle 
%             varargin  ... additional parameters for FHandle

%     Output: L         ... right hand side

function L = rhs2d(vertices,FHandle,varargin)
  
  % initialize constants
  n = size(vertices,1);
  
  % preallocate memory
  L = zeros(n,n);
  Lloc = zeros(2,2);
  
  % compute Gauss points
  [xg,wg] = gauleg(3);
  e = ones(length(wg),1);
  
  % precompute shape functions
  N = shap(xg);
  
  % assemble element contributions
  vidxi = [1 2];
  for i = 1:n-1
      vidxj = [1 2];
      for j = 1:n-1

          % compute element mapping
          
          a = vertices(vidxi(1));
          hi = vertices(vidxi(2))-a;
          xi = a + (xg+1)*hi/2;
          xi = kron(xi,e);
          a = vertices(vidxj(1));
          hj = vertices(vidxj(2))-a;
          xj = a + (xg+1)*hj/2;
          xj = kron(e,xj);
          w = kron(wg,wg);
          
          % compute load data
          FVal = FHandle(xj,xi,varargin{:});
          
          % compute element load vector
          Lloc(1,1) = sum(w.*FVal.*N(:,1))*hi*hj/4;
          Lloc(2,1) = sum(w.*FVal.*N(:,2))*hi*hj/4; 
          Lloc(1,2) = sum(w.*FVal.*N(:,3))*hi*hj/4;
          Lloc(2,2) = sum(w.*FVal.*N(:,4))*hi*hj/4; 
          
          % add contributions to global load vector
          L(vidxi,vidxj) = L(vidxi,vidxj) + Lloc;
        
          % update current element 
          vidxj = vidxj+1;
      end
      vidxi = vidxi+1;
  end

  
  % degree of freedoms
  L = L(:);
  ndof = [0;ones(n-2,1);0];  
  ndof = logical(kron(ndof,ndof));
  L = L(ndof);
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SHAP computes the values of the two-dimensional 
%     shape functions for bilinear elements  
%
%     N = shap(x)
%     
%     Input:  x ... points
%
%     Output: N ... shape functions (dim n^2x4)

function N = shap(x)

  n = size(x,1);
  
  % preallocate memory
  N = zeros(n^2,4);
  
  % compute function values
  N(:,1) = kron(1/2*(1-x),1/2*(1-x));
  N(:,2) = kron(1/2*(1+x),1/2*(1-x));
  N(:,3) = kron(1/2*(1-x),1/2*(1+x));
  N(:,4) = kron(1/2*(1+x),1/2*(1+x));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GAULEG computes Gauss quadrature points 
%
%     [x w] = gauleg(n)
%     
%     Input:  n ... number of gauss points
%
%     Output: x ... gauss points
%             w ... weights

function [x,w] = gauleg(n)
  
  % Initalize variables
  m = floor((n+1)/2);
  x  = zeros(n,1);   
  w  = zeros(n,1);
  
  for i = 1:m
    
    % Initial guess of root (starting value)
    z = cos(pi*(i-1/4)/(n+1/2));
    
    delta = 1;
    while(delta > eps)
        
      p1 = 0;
      p2 = 1;
      
      for k = 0:(n-1)

        % Computing value of n-th Legendre polynomial at point z using the
        % recursion:
        %
        %   (j+1)*P_(j+1)(z) = (2*j+1)*z*P_(j)(z)-j*P_(j-1)(z)
          
        p3 = ((2*k+1)*z*p2-k*p1)/(k+1);
        
        % Computing value of first derivative of n-th Legendre polynomial
        % at point z using the recursion:
        %
        %   (1-z^2)*P'_(j)(z) = j*[z*P_(j)(z)-P_(j-1)(z)]
        
        dp = n*(z*p3-p2)/(z^2-1);
        p1 = p2;
        p2 = p3;
        
      end    
    
      % Performing Newton update
      
      z_old = z;
      z = z_old-p3/dp;
      
      delta = abs(z-z_old);
      
    end
    
    % Compute Gauss points in [-1 1]
    x(i) = -z;
    x(n+1-i) = z;
    w(i) = 2/((1-z^2)*dp^2);
    w(n+1-i) = w(i);
    
  end
  
 return
