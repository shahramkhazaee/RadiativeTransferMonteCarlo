function obs = radiativeTransfer( physics, source, material, observation )

% discretization in packets of particles  (for optimal vectorization)
Npk = 5e4;                             % size of packets (5e4 seems optimal on my computer)
Np = ceil(source.numberParticles/Npk); % number of packets

% OBSERVATION STRUCTURE
% d         : dimension of the problem
% acoustics : true=acoustics, false=elastics
% t         : time instants
% Nt        : number of time instants
% psi       : propagation directions
% binPsi    : bins for histograms in direction
% Npsi      : number of directions
% binX      : bins for histograms in positions
% x         : sensor positions
% Nx        : number of positions
% energy    : matrix of observations size [Nx Nth Nt 2]
% dV        : small volume of domain
% dE        : energy of a single particle
obs = initializeObservation( physics, observation, Np*Npk );
energy = obs.energy;
material = prepareSigma(material);        % prepare scattering cross sections 

% loop on packages of particles
parfor ip = 1:Np

    % PARTICLES
    % N            : number of particles
    % d            : dimension of the problem
    % acoustics    : true=acoustics, false=elastics
    % x            : cartesian coordinates
    % dir          : direction of propagation
    % perp         : orthogonal to direction of propagation
    % p            : polarization (used only in elasticity)
    % meanFreePath : mean free path
    % v            : propagation velocity
    % t            : current time for the particle
    P = initializeParticle( Npk, physics, source, material );
    obsi = observeTime(obs,P,1);

    % loop on time
    for i1 = 2:obs.Nt

        % propagate particles
        P = propagateParticle(material,P,obs.t(i1));

        % observe energies (as a function of [Psi x t p])
        obsi = observeTime(obsi,P,i1);

    % end of loop on time
    end

    energy = energy + obsi.energy;

% end of loop on packages
end

% accumulate energy over all packages
obs.energy = energy;

% energy density as a function of [x t p] and [t p]
obs.energyDensity = squeeze(tensorprod(obs.dpsi',obs.energy,1));
obs.energyDomain = squeeze(tensorprod(obs.dx',obs.energyDensity,1));

% TODO
% --- check normalizations are correct by running the same simulation for
% different numbers of particles, same order of magnitude of result should
% be obtained
% --- rather than computing histograms with uniform angles, compute
% histograms with uniform angle*cosine
