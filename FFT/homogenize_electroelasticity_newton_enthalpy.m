function response = homogenize_electroelasticity_newton_enthalpy(...
    strain_macro, D_macro, geometry, material, options, GreenCurlFree, GreenDivFree)  %% FIXED-POINT SCHEME
% outputs = [strain, D, stress, E]

tic;
max_iter = options.max_iter;    % Maximum number of iterations allowed in our fixed-point scheme
TOL = options.TOL;              % Tolerance for deciding to stop the iterations
max_iter_CG = options.max_iter_CG;
TOL_CG = options.TOL_CG;
load_case_number = options.load_case_number;
n = geometry.n;
phase = geometry.phase;
h = geometry.h;
volume = geometry.volume;
id = geometry.id;
ndim = length(n);
n_phases = length(phase);

S_micro_phase = material.S_micro_phase;

n_nodes = prod(n);
T_mm = zeros(n_nodes, ndim, ndim, ndim, ndim);
T_me = zeros(n_nodes, ndim, ndim, ndim);
T_em = zeros(n_nodes, ndim, ndim, ndim);
T_ee = zeros(n_nodes, ndim, ndim);
for p = 1:n_phases
    % Derivative of stress w.r.t strain
    T_mm(phase{p},1,1,1,1) = S_micro_phase(1,1,p);
    T_mm(phase{p},1,1,2,2) = S_micro_phase(1,2,p);
    T_mm(phase{p},1,1,3,3) = S_micro_phase(1,3,p);

    T_mm(phase{p},2,2,1,1) = S_micro_phase(2,1,p);
    T_mm(phase{p},2,2,2,2) = S_micro_phase(2,2,p);
    T_mm(phase{p},2,2,3,3) = S_micro_phase(2,3,p);

    T_mm(phase{p},3,3,1,1) = S_micro_phase(3,1,p);
    T_mm(phase{p},3,3,2,2) = S_micro_phase(3,2,p);
    T_mm(phase{p},3,3,3,3) = S_micro_phase(3,3,p);

    T_mm(phase{p},2,3,2,3) = S_micro_phase(4,4,p) / 2;
    T_mm(phase{p},2,3,3,2) = S_micro_phase(4,4,p) / 2;
    T_mm(phase{p},3,2,2,3) = S_micro_phase(4,4,p) / 2;
    T_mm(phase{p},3,2,3,2) = S_micro_phase(4,4,p) / 2;

    T_mm(phase{p},1,3,1,3) = S_micro_phase(5,5,p) / 2;
    T_mm(phase{p},1,3,3,1) = S_micro_phase(5,5,p) / 2;
    T_mm(phase{p},3,1,1,3) = S_micro_phase(5,5,p) / 2;
    T_mm(phase{p},3,1,3,1) = S_micro_phase(5,5,p) / 2;

    T_mm(phase{p},1,2,1,2) = S_micro_phase(6,6,p) / 2;
    T_mm(phase{p},1,2,2,1) = S_micro_phase(6,6,p) / 2;
    T_mm(phase{p},2,1,1,2) = S_micro_phase(6,6,p) / 2;
    T_mm(phase{p},2,1,2,1) = S_micro_phase(6,6,p) / 2;

    % Derivative of stress w.r.t. electric field E
    T_me(phase{p}, 1,1, 3) = S_micro_phase(1,9,p);
    T_me(phase{p}, 2,2, 3) = S_micro_phase(2,9,p);
    T_me(phase{p}, 3,3, 3) = S_micro_phase(3,9,p);

    T_me(phase{p}, 2,3, 2) = S_micro_phase(4,8,p) / sqrt(2);
    T_me(phase{p}, 3,2, 2) = S_micro_phase(4,8,p) / sqrt(2);
    T_me(phase{p}, 1,3, 1) = S_micro_phase(5,7,p) / sqrt(2);
    T_me(phase{p}, 3,1, 1) = S_micro_phase(5,7,p) / sqrt(2);

    % Derivative of electric induction D w.r.t. strain
    T_em(phase{p}, 2, 2,3) = S_micro_phase(8,4,p) / sqrt(2);
    T_em(phase{p}, 2, 3,2) = S_micro_phase(8,4,p) / sqrt(2);
    T_em(phase{p}, 1, 1,3) = S_micro_phase(7,5,p) / sqrt(2);
    T_em(phase{p}, 1, 3,1) = S_micro_phase(7,5,p) / sqrt(2);

    T_em(phase{p}, 3, 1,1) = S_micro_phase(9,1,p);
    T_em(phase{p}, 3, 2,2) = S_micro_phase(9,2,p);
    T_em(phase{p}, 3, 3,3) = S_micro_phase(9,3,p);

    % Deriative of electric induction D w.r.t. electric field E
    T_ee(phase{p}, 1,1) = S_micro_phase(7,7,p);
    T_ee(phase{p}, 2,2) = S_micro_phase(8,8,p);
    T_ee(phase{p}, 3,3) = S_micro_phase(9,9,p);
end
T_mm = reshape(T_mm, [n, ndim, ndim, ndim, ndim]);
T_me = reshape(T_me, [n, ndim, ndim, ndim]);
T_em = reshape(T_em, [n, ndim, ndim, ndim]);
T_ee = reshape(T_ee, [n, ndim, ndim]);

strain = permute(repmat(strain_macro, [1,1, n_nodes]), [3, 1, 2]);
D = permute(repmat(D_macro, 1, n_nodes), [2, 1]);

IND_1 = 1 : (n_nodes*ndim*ndim);
IND_2 = n_nodes*ndim*ndim + (1:n_nodes *ndim);

stress = zeros(n_nodes, ndim, ndim);
E = zeros(n_nodes, ndim);
stress_fourier = zeros([n, ndim, ndim]);
E_fourier = zeros([n, ndim]);

for iter = 0 : max_iter
    % Compute dual fields according to linear constitutive law
    for p = 1:n_phases
        [stress(phase{p},:,:), E(phase{p},:)] = constitutive_law(...
            strain(phase{p},:,:), D(phase{p},:), S_micro_phase(:,:,p));
    end
    stress = reshape(stress, [n, ndim, ndim]);
    E = reshape(E, [n, ndim]);
    for i = 1:ndim
        for j = 1:ndim
            stress_fourier(:,:,:,i,j) = fftn(stress(:,:,:,i,j));
        end
        E_fourier(:,:,:,i) = fftn(E(:,:,:,i));
    end
    stress = reshape(stress, [n_nodes, ndim, ndim]);
    E = reshape(E, [n_nodes, ndim]);

    lhs = @(dgrad) stiff_func3D(dgrad, IND_1, IND_2, T_mm, T_me, T_em, T_ee, GreenCurlFree, GreenDivFree);            % LEFT-HAND SIDE
    rhs = residual_func3D(stress_fourier, E_fourier, GreenCurlFree, GreenDivFree);                 % RIGHT-HAND SIDE
    residual = norm(rhs, 2) / sqrt(length(rhs));
%     fprintf('Step %d: residual = %e, ', iter, residual);

    if residual < TOL
        fprintf('%06d Load%2d- Converged at iteration = %d \n', id, load_case_number, iter);
        break;
    end
    dgrad = pcg(lhs, rhs, TOL_CG, max_iter_CG);
    % dgrad = bicgstabl(lhs, rhs, TOL_CG, max_iter_CG);

    strain = strain + reshape(dgrad(IND_1), [n_nodes, ndim, ndim]);
    D = D + reshape(dgrad(IND_2), [n_nodes, ndim]);
end

if residual > TOL
    fprintf('%06d Load%2d- Not Converged after %d iterations! \n', id, load_case_number, iter);
    response = 'failed';
    return;
end

stress = reshape(stress, [n, ndim, ndim]);
E = reshape(E, [n, ndim]);

% stress_macro{q} = h * squeeze(sum(stress, 1:3)) / volume;
% E_macro{q} = h * squeeze(sum(E, 1:3)) / volume;
response.id = id;
response.stress_macro = h * squeeze(sum(stress, 1:3)) / volume;
response.E_macro = h * squeeze(sum(E, 1:3)) / volume;
response.time = toc;
end

%% HELPER FUNCTIONS
function [stress, E] = constitutive_law(strain, D, S)
stress = zeros(size(strain));
E = zeros(size(D));

stress(:, 1,1) = S(1,1) * strain(:, 1,1) + S(1,2) * strain(:, 2,2) + S(1,3) * strain(:, 3,3) ...
    + S(1,9) * D(:, 3);
stress(:, 2,2) = S(2,1) * strain(:, 1,1) + S(2,2) * strain(:, 2,2) + S(2,3) * strain(:, 3,3) ...
    + S(2,9) * D(:, 3);
stress(:, 3,3) = S(3,1) * strain(:, 1,1) + S(3,2) * strain(:, 2,2) + S(3,3) * strain(:, 3,3) ...
    + S(3,9) * D(:, 3);

stress(:, 2,3) = S(4,4) * strain(:, 2,3) + S(4,8) / sqrt(2) * D(:, 2) ;
stress(:, 1,3) = S(5,5) * strain(:, 1,3) + S(5,7) / sqrt(2) * D(:, 1);
stress(:, 1,2) = S(6,6) * strain(:, 1,2);

stress(:, 3,2) = stress(:, 2,3);
stress(:, 3,1) = stress(:, 1,3);
stress(:, 2,1) = stress(:, 1,2);

E(:, 1) = S(7,5) * sqrt(2) * strain(:, 1,3) + S(7,7) * D(:, 1);
E(:, 2) = S(8,4) * sqrt(2) * strain(:, 2,3) + S(8,8) * D(:, 2);
E(:, 3) = S(9,1) * strain(:, 1,1) + S(9,2) * strain(:, 2,2) + S(9,3) * strain(:, 3,3) ...
    + S(9,9) * D(:, 3);
end

function v = residual_func3D(stress_fourier, E_fourier, GreenCurlFree, GreenDivFree)
n = zeros(1,3);
[n(1), n(2), n(3), ~] = size(GreenCurlFree);
ndim = 3;
v1 = zeros([n, ndim, ndim]);
v2 = zeros([n, ndim]);

R = zeros([n, ndim, ndim]);
for i = 1:ndim
    for j = 1:ndim
        R(:,:,:,i,j) = squeeze(sum( GreenCurlFree(:,:,:,j,:) .* stress_fourier(:,:,:,i,:), ndim+2));
    end
end

for i = 1:ndim
    for j = 1:ndim
        v1(:,:,:,i,j) = -0.5 * ifftn( R(:,:,:,i,j) + R(:,:,:,j,i) );
    end
    v2(:,:,:,i) = -ifftn( squeeze( sum( squeeze(GreenDivFree(:,:,:,i,:)) .* E_fourier, 4) ) );
end

v = vertcat(v1(:), v2(:));
end

function v = stiff_func3D(dgrad, IND_1, IND_2, T_mm, T_me, T_em, T_ee, GreenCurlFree, GreenDivFree)
n = zeros(1,3);
[n(1), n(2), n(3), ~] = size(GreenCurlFree);
ndim = 3;
dstrain = reshape(dgrad(IND_1), [n, ndim, ndim]);
dD = reshape(dgrad(IND_2), [n, ndim]);

dT1 = zeros([n, ndim, ndim]);
dT2 = zeros([n, ndim]);
for i = 1:ndim
    for j = 1:ndim
        dT1(:,:,:,i,j) = fftn( squeeze(sum( squeeze(T_mm(:,:,:,i,j,:,:)) .* dstrain, [4, 5] ) ) ...
            + squeeze(sum( squeeze(T_me(:,:,:,i,j,:)) .* dD, 4)) );
    end
    dT2(:,:,:,i) = fftn( squeeze( sum( squeeze(T_em(:,:,:,i,:,:)) .* dstrain, [4, 5]) ) ...
        + squeeze( sum( squeeze(T_ee(:,:,:,i,:)) .* dD, 4) ) );
end
v1 = zeros([n, ndim, ndim]);
v2 = zeros([n, ndim]);
R = zeros([n, ndim, ndim]);
for i = 1:ndim
    for j = 1:ndim
        R(:,:,:,i,j) = squeeze( sum(GreenCurlFree(:,:,:,j,:) .* dT1(:,:,:,i,:), ndim+2) );
    end
end
for i = 1:ndim
    for j = 1:ndim
        v1(:,:,:,i,j) = 0.5 * ifftn( R(:,:,:,i,j) + R(:,:,:,j,i) );
    end
    v2(:,:,:,i) = ifftn( squeeze(sum(squeeze(GreenDivFree(:,:,:,i,:)) .* dT2, 4)) );
end

v = vertcat(v1(:), v2(:));
end