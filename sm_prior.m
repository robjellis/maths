function [xadj xhat] = sm_prior(x,param)

% smoothness priors function
%
% code from: http://www.mathworks.com/matlabcentral/newsreader/view_thread/22445
% by ghassan@s2.chalmers.se ("Ghassan Hamarneh")
%

% Smoothness priors¹ based method, it is easy to control with one
% parameter and it can be expressed using matrix notation, thus it's
% easy to implement with matlab. If we base the smoothing to second
% difference matrix, the magic matlab lines will be:
%
% ¹W. Gersch, Smoothness priors , in New Directions in Time Series
% Analysis, Part II, pp. 113 146. Springer-Verlag, 1991.

x = x(:);

N=numel(x); % Number of points

e=ones(N,1); % building of the 2nd diff. matrix
L=spdiags([e -2*e e], 0:2, N-2, N);

alpha=param; % the control parameter, larger -> smoother

%xhat=inv(speye(N,N)+alpha^2*(L'*L))*x; % Estimation/smoothing
% inv is slow. Instead of inv(A)*b, use A\b (per MATLAB)

xhat= (speye(N,N)+alpha^2*(L'*L)) \ x; % Estimation/smoothing

xadj = x - xhat; % the adjusted x-values
figure(100)
plot(1:N,x,'b',1:N,xhat,'k') % Visualize results
