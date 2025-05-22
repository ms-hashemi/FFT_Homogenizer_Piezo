function response = homogenize_electroelasticity_newton_energy_sym_efficient(...
    strain_macro, D_macro, geometry, material, options, GreenCurlFree, GreenDivFree)  %% FIXED-POINT SCHEME

tic;
max_iter = options.max_iter;    % Maximum number of iterations allowed in our fixed-point scheme
TOL = options.TOL;              % Tolerance for deciding to stop the iterations
max_iter_CG = options.max_iter_CG;
TOL_CG = options.TOL_CG;
n = geometry.n;
phase = geometry.phase;
h = geometry.h;
volume = geometry.volume;
id = geometry.id;
ndim = length(n);
sdim = ndim * (ndim + 1)/2;
n_phases = length(phase);

S_micro_phase = material.S_micro_phase;

n_nodes = prod(n);
strain = permute(repmat(strain_macro, [1, n_nodes]), [2, 1]);
D = permute(repmat(D_macro, 1, n_nodes), [2, 1]);

IND_mech = 1 : (n_nodes * sdim);
IND_elec = n_nodes * sdim + (1:n_nodes * ndim);

stress = zeros(size(strain));
E = zeros(n_nodes, ndim);
stress_fourier = zeros([n, sdim]);
E_fourier = zeros([n, ndim]);

lhs = @(dgrad) linear_operator(dgrad, IND_mech, IND_elec, ...
    GreenCurlFree, GreenDivFree, phase, S_micro_phase);            % LEFT-HAND SIDE
for iter = 1 : max_iter
    % Compute dual fields according to linear constitutive law
    for p = 1:n_phases
        [stress(phase{p},:), E(phase{p},:)] = constitutive_law(...
            strain(phase{p},:), D(phase{p},:), S_micro_phase(:,:,p));
    end
    
    stress = reshape(stress, [n, sdim]);
    E = reshape(E, [n, ndim]);

    for q = 1:sdim
        stress_fourier(:,:,:,q) = fftn(stress(:,:,:,q));
    end
    for i = 1:ndim
        E_fourier(:,:,:,i) = fftn(E(:,:,:,i));
    end
    stress = reshape(stress, [n_nodes, sdim]);
    E = reshape(E, [n_nodes, ndim]);

    rhs = residual_func(stress_fourier, E_fourier, GreenCurlFree, GreenDivFree);    % RIGHT-HAND SIDE
    residual = norm(rhs, 2) / sqrt(length(rhs));
%     fprintf('Step %d: Residual = %e ', iter, residual);

    if residual < TOL
        fprintf('\nConvergenced reached at iteration step = %d! \n', iter);
        break;              
    end
    dgrad = pcg(lhs, rhs, TOL_CG, max_iter_CG);

    strain = strain + reshape(dgrad(IND_mech), [n_nodes, sdim]);
    D = D + reshape(dgrad(IND_elec), [n_nodes, ndim]);
end

if residual > TOL
%     fprintf('%06d Load%2d- Not Converged after %d iterations! \n', id, load_case_number, iter);
    fprintf('%s Load%2d- Not Converged after %d iterations! \n', id, load_case_number, iter);
    response = 'failed';
    return;
end

response.id = id;
response.stress_macro = h * squeeze(sum(stress, 1))' ./ volume;
response.E_macro = h * squeeze(sum(E, 1))' ./ volume;
response.time = toc;
% response.stress = stress;
% response.E = E;
% response.strain = strain;
% response.D = D;
% save('test.mat', 'stress', 'E', 'strain', 'D');
end

%% HELPER FUNCTIONS
function [stress, E] = constitutive_law(strain, D, moduli)
stress = zeros(size(strain));
E = zeros(size(D));

stress(:, 1) = moduli(1,1) * strain(:, 1) + moduli(1,2) * strain(:, 2) ...
    + moduli(1,3) * strain(:, 3) + moduli(1,9) * D(:, 3);
stress(:, 2) = moduli(2,1) * strain(:, 1) + moduli(2,2) * strain(:, 2) ...
    + moduli(2,3) * strain(:, 3) + moduli(2,9) * D(:, 3);
stress(:, 3) = moduli(3,1) * strain(:, 1) + moduli(3,2) * strain(:, 2) ...
    + moduli(3,3) * strain(:, 3) + moduli(3,9) * D(:, 3);

stress(:, 4) = moduli(4,4) * strain(:, 4) + moduli(4,8) / sqrt(2) * D(:, 2);  % component stress(2,3) --> stress(4)
stress(:, 5) = moduli(5,5) * strain(:, 5) + moduli(5,7) / sqrt(2) * D(:, 1);  %
stress(:, 6) = moduli(6,6) * strain(:, 6);

E(:, 1) = moduli(7,5) * sqrt(2) * strain(:, 5) + moduli(7,7) * D(:, 1);
E(:, 2) = moduli(8,4) * sqrt(2) * strain(:, 4) + moduli(8,8) * D(:, 2);
E(:, 3) = moduli(9,1) * strain(:, 1) + moduli(9,2) * strain(:, 2) ...
    + moduli(9,3) * strain(:, 3) + moduli(9,9) * D(:, 3);
end

function v = residual_func(stress_fourier, E_fourier, GreenCurlFree, GreenDivFree)
ndim = 3;
n = zeros(1, ndim);
[n(1), n(2), n(3), sdim] = size(GreenCurlFree);
IND = [1, 6, 5; 
    6, 2, 4; 
    5, 4, 3];
IND_sym = [1  1;  2  2;  3  3;  
    2  3;  1  3;  1  2];

R = zeros([n, ndim, ndim]);
for i = 1:ndim
    for j = 1:ndim
        R(:,:,:,i,j) = squeeze(sum( GreenCurlFree(:,:,:,IND(i,:)) .* stress_fourier(:,:,:,IND(j,:)) , ndim+1));
    end
end

v1 = zeros([n, sdim]);
for q = 1:sdim
    v1(:,:,:,q) = -ifftn( R(:,:,:,IND_sym(q,1), IND_sym(q,2)) + R(:,:,:,IND_sym(q,2), IND_sym(q,1)) );
end

v2 = zeros([n, ndim]);
for i = 1:ndim
    v2(:,:,:,i) = -ifftn( squeeze( sum( GreenDivFree(:,:,:, IND(i,:)) .* E_fourier, ndim+1) ) );
end

v = vertcat(v1(:), v2(:));
end

function v = linear_operator(dgrad, IND_mech, IND_elec, ...
    GreenCurlFree, GreenDivFree, phase, moduli)

ndim = 3;
n = zeros(1,ndim);
[n(1), n(2), n(3), sdim] = size(GreenCurlFree);
n_nodes = prod(n);

dstrain = reshape(dgrad(IND_mech), [n_nodes, sdim]);
dD = reshape(dgrad(IND_elec), [n_nodes, ndim]);

dT_mech = zeros([n_nodes, sdim]);
dT_elec = zeros([n_nodes, ndim]);

n_phases = length(phase);
for p = 1:n_phases
    [dT_mech(phase{p}, :), dT_elec(phase{p},:)] ...
        = constitutive_law(dstrain(phase{p}, :), dD(phase{p},:), moduli(:,:,p));
end
dT_mech = reshape(dT_mech, [n, sdim]);
dT_elec = reshape(dT_elec, [n, ndim]);

for q = 1:sdim
    dT_mech(:,:,:,q) = fftn( dT_mech(:,:,:,q) );
end
for i = 1:ndim
    dT_elec(:,:,:,i) = fftn( dT_elec(:,:,:,i) );
end

IND = [1,  6,  5;
       6,  2,  4;
       5,  4,  3];
IND_sym = [1  1; 2  2; 3  3;
    2  3; 1  3; 1  2];

R = zeros([n, ndim, ndim]);
for i = 1:ndim
    for j = 1:ndim
        R(:,:,:,i,j) = squeeze( sum( GreenCurlFree(:,:,:,IND(i,:)) .* dT_mech(:,:,:,IND(j,:)), ndim+1) );
    end
end

v1 = zeros([n, sdim]);
for q = 1:sdim
    v1(:,:,:,q) = ifftn( R(:,:,:, IND_sym(q,1), IND_sym(q,2)) + R(:,:,:, IND_sym(q,2), IND_sym(q,1)) );
end
v2 = zeros([n, ndim]);
for i = 1:ndim
    v2(:,:,:,i) = ifftn( squeeze( sum( GreenDivFree(:,:,:,IND(i,:)) .* dT_elec, ndim+1 ) ) );
end

v = vertcat(v1(:), v2(:));
end