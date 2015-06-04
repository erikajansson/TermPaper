function [] = mySlope(axes_handle,loc,pos,rate,cstr)

% ADD_SLOPE Add slope triangle to current plot.
%
%   ADD_SLOPE(AXES_HANDLE,LOC,POS,RATE) adds a slope triangle with convergence
%   rate P at the location LOC with postion POS to the plot specified 
%   by AXES_HANDLE.
%

% initialize constants
F_SIZE = 12;  % Font size

% get axes limits
Xlim = get(axes_handle,'XLim');

% Compute triangle location
if (pos==1)      
    
    % compute vertices of triangle
    xr = loc(1);
    xl = loc(1)*exp(log(Xlim(2)/Xlim(1))*0.2);
    xmid = sqrt(xl*xr);
    yr = loc(2);
    yl = yr*(xl/xr)^rate;
    valign = 'bottom';
    %valign = 'top';
    ymid = yl;
    
else 
       
    % compute vertices of triangle  
    xr = loc(1)*exp(log(Xlim(2)/Xlim(1))*0.2);
    xl = loc(1);
    xmid = sqrt(xl*xr);
    yl = loc(2);
    yr = yl*(xr/xl)^rate;
    valign = 'top';
    %valign = 'bottom';
    ymid = yl;
    
end
  
  
% add slope triangle to current plot
h = plot([xl xr xr xl],[yl yl yr yl],cstr);
cstr = get(h,'Color'); 
text(xmid,ymid,sprintf('s = %1.1f',rate), ...
       'FontSize',F_SIZE, ...
       'HorizontalAlignment','center', ...
       'VerticalAlignment',valign, ...
       'Color',cstr);
   
   
  
  
return



