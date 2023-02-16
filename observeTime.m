function obs = observeTime(obs,P,it)

% compute energy normalization
dE = obs.dE./(obs.dpsi'*obs.dx);

% compute radius and angle between direction of propagation and position
% check normalization of dir
r = sqrt(sum(P.x.^2,2));
cospsi = dot(P.x,P.dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
psi = acos(cospsi);

% accumulate energies
ind = P.p;
Ep = histcounts2( psi(ind), r(ind), obs.binPsi, obs.binX ).*dE;
Es = histcounts2( psi(~ind), r(~ind), obs.binPsi, obs.binX ).*dE;
if it==1
    obs.energy(:,:,1) = Ep;
    obs.energy(:,:,2) = Es;
else
    obs.energy(:,:,it,1) = obs.energy(:,:,it,1) + Ep;
    obs.energy(:,:,it,2) = obs.energy(:,:,it,2) + Es;
end