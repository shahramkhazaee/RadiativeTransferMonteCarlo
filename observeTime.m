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
E = histcounts2( psi, r, obs.binPsi, obs.binX ).*dE;
if it==1
    obs.energy(:,:,1) = E;
else
    obs.energy(:,:,it) = obs.energy(:,:,it) + E;
end