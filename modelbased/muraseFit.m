function [par,fit,cnum] = muraseFit(time,plasma,tissue,num_par)

% default to fitting with k1,k2,Vp
if (nargin < 4) num_par = 3; end;

if (num_par < 2 | num_par > 3) error('muraseFit :: number of parameters must be 2 or 3'); end;

if (min(size(time)) ~= 1) error('muraseFit :: time and plasma curves must be vectors'); end;

if (size(time,1) == 1) time = time'; end;
if (size(plasma,1) == 1) plasma = plasma'; end;
if (size(tissue,1) == 1) tissue = tissue'; end;

if (size(time) ~= size(plasma)) error('muraseFit :: time and plasma curves are inconsistent'); end;
if (size(time) ~= size(tissue)) error('muraseFit :: time and tissue curves are inconsistent'); end;

A = zeros([length(time) num_par]);

if (num_par == 2)
    A(:,1) = cumtrapz(time,plasma);
    A(:,2) = -cumtrapz(time,tissue);
elseif (num_par == 3)
    A(:,1) = cumtrapz(time,plasma);
    A(:,2) = -cumtrapz(time,tissue);
    A(:,3) = plasma;  
end;

% solve design matrix via SVD
% close all
warning off all

par = A\tissue;

fit = A*par;

% plot(tissue,'x')
% hold on
% plot(fit)        
% pause

if (num_par == 3)
    % compute k1
    par(1) = par(1) - par(2)*par(3);
end;

par = par';

% output design matrix condition number
if (nargout == 3) cnum = cond(A); end;

return;


