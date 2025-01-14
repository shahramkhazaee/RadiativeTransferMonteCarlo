function [E,E_diff] = Energy_Paasschens(d,a,b)
% Function calculating the solution of the RTE of scalar waves propagating 
% in an isotropic scattering random medium excited by an isotropic (explosion) 
% source with the Diffusion approx. in 2D & 3D based on Paaschens's 1997 paper. 
% Note: the solution is exact for 2D and approximate for 3D.

% inputs (normalized parameters)
% a (normalized time)      : Sigma*v*t
% b (normalized distance)  : r*Sigma, i.e. source-station distance * scattering cross-section
% d (dimension of domain)  : either 2 or 3

% outputs
% E : normalized energy density 
    % in 2D : E : r*E/Sigma
    % in 3D : E : r^2*E/Sigma
% E_diff : diffusion approximation (same normalization used here)

[m,n] = size(a);
if m>n, a=a'; end
[m,n] = size(b);
if m<n, b=b'; end

if d==2
    E = b./(2*pi*sqrt(a.^2-b.^2)).*exp(sqrt(a.^2-b.^2)-a).*heaviside(a-b);
    
    aa = repmat(a,[length(b) 1]);
    ind = (a==b);
    if ~isempty(E(ind))
        E(ind) = E(ind) + exp(-aa(ind))/(2*pi);
    end

    if nargout==2
        E_diff = b*d.*exp(-d*b.^2./(4*a))./(4*pi*a);
    end
elseif d==3
    fun = @(x) exp(x).*sqrt(1+2.026./x);
    E = b.^2.*(1-(b./a).^2).^(1/8).*exp(-a).*fun(a.*(1-(b./a).^2).^(3/4)).*heaviside(a-b)./(4*pi*a/3).^(3/2);
    
    aa = repmat(a,[length(b) 1]); bb = repmat(b,[1 length(a)]);
    ind = (a==b);
    if ~isempty(E(ind))
        E(ind) = E(ind) + exp(-aa(ind))./(4*pi*bb(ind));
    end

    if nargout==2
        E_diff = b.^2.*exp(-d*b.^2./(4*a))./(4*pi*a/d).^(3/2);
    end
else
    disp('Dimension "d" should be either 2 or 3 !')
end