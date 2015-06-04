%  BS_EXCHANGE computes value of a exchange option 
%     with the Black-Scholes model using
%     analytic formulas 
%    
%     [P] = bs_exhange([S1 S2],T,[sigma1 sigma2],rho)
%     
%     Input:  S1,S2  ... stock prices at time 0
%             T      ... maturity
%             sigma1 ... volatility of the first GBM
%             sigma2 ... volatility of the second GBM
%             rho    ... correlation 
%             a      ... payoff constant
%             b      ... payoff constant
%
%     Output: P     ... option price at time 0

function P = bs_exchange(S,T,sigma,rho,a,b)

S1 = a*S(:,1);
S2 = b*S(:,2);
s = S1./S2;

sgm = sqrt(sigma(1)^2+sigma(2)^2-2*sigma(1)*sigma(2)*rho);
d1 = (log(S1./S2)+sgm^2*T/2)/(sgm*sqrt(T));
d2 = d1-sgm*sqrt(T);

P = S1.*normcdf(d1) - S2.*normcdf(d2);

return


