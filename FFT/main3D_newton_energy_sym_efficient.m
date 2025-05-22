%% PARALLELIZATION

parallel_size = 9;
delete(gcp('nocreate'))
p = parpool(parallel_size);

current_path = cd;
current_path = string(current_path);
% path_input_images = fullfile('..', '..', '..', '..', 'ISU', 'Com S 527 Concurrent programming', 'Project', 'ML', 'data', 'deg', 'Selected');
% path_input_images = fullfile('..', 'deg', 'Selected');
path_input_images = fullfile('.', 'test');
% microstructure_ids = importdata('job_list.txt');
file = fopen('job_list.txt', 'r');
number = 0;
while ~feof(file)
    str = fgetl(file);
    number = number + 1;
end
microstructure_ids = cell(number, 1);
frewind(file);
number = 1;
while ~feof(file)
    microstructure_ids{number} = fgetl(file);
    number = number + 1;
end
microstructure_ids_desired = 1;
parallel_future_objects(1:9*length(microstructure_ids_desired)) = parallel.FevalFuture;
object_counter = 1;

mkdir('Results_test');
log_file_id = fopen('log_test.txt', 'a');

%% REGULAR MESH FOR 3D DOMAIN
ndim = 3;             % Problem dimension
L = [1, 1, 1/5];        % RVE occupies the domain [-L(1), L(1)] x [-L(2), L(2)] x [-L(3), L(3)]
volume = prod(2*L);

% if ~exist('n', 'var')
%     n = repmat(64, [1, 3]);
% elseif length(n) == 1
%     n = repmat(n,   [1, 3]);
% end
% % [phase_IND, n] = read_microimage_resize();
% h = prod(2 * L./n);     % size of one cell in the FFT method
% scale = L./pi;          % for scaling the wave frequency

% Load microstructure data from a previously saved MATLAB data file
% micro_size = [128, 128, 64];
% file_name = ['microstructure_image_', sprintf('%dx%dx%d', micro_size(1), micro_size(2), micro_size(3)), '.mat'];
% load(file_name);
% n = size(phase_IND);

for id = microstructure_ids_desired
    % Get microstructure data from its sliced images
    % [phase_IND, n] = read_microimage_resize('000009', 'slice_', 0, 'png', [1, 150; 1, 150; 1, 75],[64, 64, 64]);
    % [phase_IND, n] = read_microimage_resize('000010.zip', 'slice_', 0, 'png', [1, 150; 1, 150; 1, 75], [64, 64, 64]);

%     [phase_IND, n] = read_microimage_resize(fullfile(path_input_images, microstructure_ids(id)), 'slice_', 0, 'png', [1, 75; 1, 75; 1, 36], [64, 64, 64]);
%     [phase_IND, n] = read_microimage_resize([fullfile(path_input_images, microstructure_ids{id}) '.zip'], 'slice_', 3, 'bmp', [16, 165; 16, 165; 16, 165], [64, 64, 64]);
    [phase_IND, n] = read_microimage_resize([fullfile(path_input_images, microstructure_ids{id})], 'slice_', 3, 'png', [16, 165; 16, 165; 16, 165], [150, 150, 150]);
    if(phase_IND == -1)
        fprintf(log_file_id, '%s- Bad input!\n', microstructure_ids{id});
%         break
        continue
    end
    h = prod(2 * L./n);     % size of one cell in the FFT method
    scale = L./pi;          % for scaling the wave frequency
    
    % Generate the mesh of grid points and wavenumbers.
    x1D = cell(1,ndim);     % one-dimensional coordinates in each direction
    xi1D = cell(1,ndim);    % one-dimensional Fourier frequencies in each direction
    for a = 1:ndim
        if rem(n(a),2) == 1
            x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2));
            xi1D{a} = [0 : fix(n(a)/2), -fix(n(a)/2) : -1] / scale(a);
        else
            x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2)-1);
            xi1D{a} = [0 : fix(n(a)/2)-1, 0, -fix(n(a)/2)+1 : -1] / scale(a);
        end
    end

    x = zeros([n, ndim]);       % three-dimensional coordinates (x)
    [x(:,:,:,1), x(:,:,:,2),x(:,:,:,3)] = ndgrid(x1D{1}, x1D{2}, x1D{3});
    xi = zeros([n, ndim]);      % three-dimensional Fourier frequencies (xi)
    [xi(:,:,:,1), xi(:,:,:,2),xi(:,:,:,3)] = ndgrid(xi1D{1}, xi1D{2}, xi1D{3});

    % % FIBER INCLUSION
    % phase_IND = ones(n);                                % matrix phase is indicated by integer 1, IND stands for indices
    % volume_frac = 0.6;
    % radius = sqrt(4 * L(1) * L(2) * volume_frac / pi);  % Fiber volume = pi*r^2 * L(3), RVE volume = 8*prod(L)
    % phase_IND(sum(x(:,:,:,1:2).^2, ndim+1) <= radius * radius) = 2;
    % 
    % xx = reshape(x, [], 3);
    % figure
    % scatter3(xx(:,1), xx(:,2), xx(:,3), 2, phase_IND(:));

    %% DEFINE GREEN OPERATOR IN FOURIER SPACE
    %===========================================================
    % Constants in Brenner's paper
    % C11_phase = [8.0,  154.837];        % lambda + 2*mu
    % C12_phase = [4.4,   83.237];        % lambda
    % C13_phase = [4.4,   82.712];        % lambda
    % C33_phase = [8.0,  131.39];         % lambda + 2*mu
    % C44_phase = 2 * [1.8,  25.696];         % 2*mu
    % C66_phase = 2 * [1.8,  35.800];         % 2*mu
    % 
    % e31_phase = [0,  -2.120582];
    % e33_phase = [0,   9.521830];
    % e15_phase = sqrt(2) * [0,   9.34959];
    % 
    % gamma11_phase = [0.0372,   4.065];
    % gamma33_phase = [0.0372,   2.079];
    %===========================================================

    %===========================================================
    % Our material properties, PDMS and Poled BTO
    C11_phase = [ 2.2839, 158.0];        % lambda + 2*mu
    C12_phase = [ 2.1580,  69.1];        % lambda
    C13_phase = [ 2.1580,  67.5];        % lambda
    C33_phase = [ 2.2839, 150.0];         % lambda + 2*mu
    C44_phase = 2 * [ 0.0629,  45.1];     % 2*mu
    C66_phase = 2 * [ 0.0629,  44.6];     % 2*mu
    % Our material properties, PDMS and Poled BTO
    C11_phase = [ 2.326, 159.0]; % lambda + 2*mu
    C12_phase = [ 2.200,  69.0]; % lambda
    C13_phase = [ 2.200,  69.0]; % lambda
    C33_phase = [ 2.3260, 159.0]; % lambda + 2*mu
    C44_phase = 2 * [ 0.063, 45.0]; % 2*mu
    C66_phase = 2 * [ 0.063, 45.0]; % 2*mu
    
    e31_phase = [ 0,  -3.14];
    e33_phase = [ 0,   14.5];
    e15_phase = [ 0,  10.9];

    epsilon_r_11_phase = [2.5, 1000]; % Relative permittivities; gamma = epsilon_r*epsilon_0 = epsilon_r*8.854e-3
    epsilon_r_33_phase = [2.5,  910]; % Relative permittivities; gamma = epsilon_r*epsilon_0 = epsilon_r*8.854e-3
    gamma11_phase = 8.854e-3.*epsilon_r_11_phase; % nC/Vm = nF/m
    gamma33_phase = 8.854e-3.*epsilon_r_33_phase; % nC/Vm = nF/m
    %===========================================================

    n_phases = 2;
    phase = cell(1,2);
    phase{1} = (phase_IND == 1);
    phase{2} = (phase_IND == 2);

    C_micro_phase = zeros(9, 9, n_phases);
    S_micro_phase = zeros(9, 9, n_phases);

    for p = 1:n_phases
        C_mm = [C11_phase(p), C12_phase(p), C13_phase(p), 0, 0, 0; ...
            C12_phase(p), C11_phase(p), C13_phase(p), 0, 0, 0; ...
            C13_phase(p), C13_phase(p), C33_phase(p), 0, 0, 0; ...
            0, 0, 0, C44_phase(p), 0, 0; ...
            0, 0, 0, 0, C44_phase(p), 0; ...
            0, 0, 0, 0, 0, C66_phase(p)];
        C_me = [0, 0, -e31_phase(p); ...
            0, 0, -e31_phase(p); ...
            0, 0, -e33_phase(p); ...
            0, -e15_phase(p), 0; ...
            -e15_phase(p), 0, 0; ...
            0, 0, 0];
        C_em = [0, 0, 0, 0, e15_phase(p), 0; ...
            0, 0, 0, e15_phase(p), 0, 0; ...
            e31_phase(p), e31_phase(p), e33_phase(p), 0, 0, 0];
        C_ee = [gamma11_phase(p), 0, 0; ...
            0, gamma11_phase(p), 0; ...
            0, 0, gamma33_phase(p)];
        C_micro_phase(:,:,p) = [C_mm, C_me; ...
                                C_em, C_ee];

        S_micro_phase(1:6,1:6, p) = C_mm - C_me * (C_ee \ C_em);
        S_micro_phase(1:6,7:9, p) = C_me / C_ee;
        S_micro_phase(7:9,1:6, p) = -C_ee \ C_em;
        S_micro_phase(7:9,7:9, p) = inv(C_ee);
    end

    material.S_micro_phase = S_micro_phase;
    optionss.max_iter = 50;
    optionss.TOL = 1e-10;
    optionss.max_iter_CG = 250;
    optionss.TOL_CG = 5e-6;

    geometry.n = n;
    geometry.phase = phase;
    geometry.h = h;
    geometry.volume = volume;
    geometry.id = microstructure_ids(id);

    %% DEFINE GREEN OPERATOR IN FOURIER SPACE
    xi_power2 = sum(xi.^2, ndim + 1);
    IND_G = [1  1;  2  2;  3  3;  
             2  3;  1  3;  1  2];
    GreenCurlFree = xi(:,:,:,IND_G(:,1)) .* xi(:,:,:,IND_G(:,2)) ./ sum(xi.^2, 4);
    nan_index = isnan(GreenCurlFree);
    inf_index = isinf(GreenCurlFree);
    GreenCurlFree(nan_index) = 0;
    GreenCurlFree(inf_index) = 0;

    I_vect = [1,  1,  1,  0,  0,  0]';
    GreenDivFree = permute(repmat(I_vect, [1, n]), [2:4, 1]) - GreenCurlFree;
    GreenDivFree(nan_index) = 0;
    GreenDivFree(inf_index) = 0;
    sdim = ndim * (ndim + 1);

    %% COMPUTE MACROSCOPIC DUAL VARIABLE ACCORDING TO UNIT LOADINGS
    
%     stress_macro = cell(1,9);
%     E_macro = cell(1,9);
%     parpool(3);
%     tic
    for q = 1:9
%         if q == 2
%             return
%         end
        strain_macro = zeros(6, 1);
        D_macro = zeros(3,1);
        switch q
            case 1
                strain_macro(1) = 1;
            case 2
                strain_macro(2) = 1;
            case 3
                strain_macro(3) = 1;
            case 4
                strain_macro(4) = 1;
            case 5
                strain_macro(5) = 1;
            case 6
                strain_macro(6) = 1;
            case 7
                D_macro(1) = 1;
            case 8
                D_macro(2) = 1;
            case 9
                D_macro(3) = 1;
        end
%         [strain, D, stress, E] = homogenize_electroelasticity_newton_energy_sym_efficient(...
%             strain_macro, D_macro, geometry, material, optionss, GreenCurlFree, GreenDivFree);
% 
%         stress_macro{q} = h * squeeze(sum(stress, 1))' ./ volume;
%         E_macro{q} = h * squeeze(sum(E, 1))' ./ volume;
        optionss.load_case_number = q;
        parallel_future_objects(object_counter) = parfeval(gcp('nocreate'), @homogenize_electroelasticity_newton_energy_sym_efficient, 1, ...
            strain_macro, D_macro, geometry, material, optionss, GreenCurlFree, GreenDivFree);
        object_counter = object_counter + 1;
    end
    if rem((object_counter - 1)/9, 8) == 0
        index = 1;
        for j = (object_counter - 1)/9-7:1:(object_counter - 1)/9
            afterall_parallel_future_objects(index) = afterAll(parallel_future_objects(9*(j-1)+1:9*(j-1)+9), @(r) calculate_save_homogenized_properties(r, log_file_id, current_path), 0);
            index = index + 1;
        end

        wait(afterall_parallel_future_objects);
    elseif length(microstructure_ids_desired) - rem(length(microstructure_ids_desired), 8) + 1 == (object_counter - 1)/9
        index = 1;
        for j = (object_counter - 1)/9:1:(object_counter - 1)/9 + rem(length(microstructure_ids_desired), 8) - 1
            afterall_parallel_future_objects2(index) = afterAll(parallel_future_objects(9*(j-1)+1:9*(j-1)+9), @(r) calculate_save_homogenized_properties(r, log_file_id, current_path), 0);
            index = index + 1;
        end

        wait(afterall_parallel_future_objects2);
    end
%     toc
end

% afterall_parallel_future_objects(1:length(microstructure_ids_desired)) = parallel.AfterAllFuture;
% for i = 1:length(microstructure_ids_desired)
%     afterall_parallel_future_objects(i) = afterAll(parallel_future_objects(9*(i-1)+1:9*(i-1)+9), @calculate_save_homogenized_properties, 0);
% end
% 
% wait(afterall_parallel_future_objects);

fclose(log_file_id);

%% Calculating the homogenized or effective properties
function calculate_save_homogenized_properties(responses, log_file_id, current_path)

for i = 1:9
    if strcmp(responses(i), 'failed')
%         fprintf(log_file_id, '\n%06d- No results due to at least one convergence failure in a load case! \n', responses(i).id);
        fprintf(log_file_id, '%s- No results due to at least one convergence failure in a load case! \n', responses(i).id{1});
        return;
    end
end

% DERIVE MACROSCOPIC LINEAR CONSTITUTIVE MATRIX (ENERGY FORM)

% In this code, we have constructed the C_macro and S_macro according to
% Brenner's paper. That is, the strain vector column is 
%
% strain-column = [ epsilon(1,1),  epsilon(2,2),  epsilon(3,3), 
%  sqrt(2) * epsilon(2,3), sqrt(2) * epsilon(1,3), sqrt(2) * epsilon(1,2) ]
%
% The stress column is then
% stress-column = [ sigma(1,1),  sigma(2,2),  sigma(3,3),
% sqrt(2) * sigma(2,3), sqrt(2) * sigma(1,3), sqrt(2) * sigma(1,2) ]

S_macro = zeros(9,9);
% IND_diag = [1, 1; 2, 2; 3 3];
% IND_offdiag = [2, 3; 1, 3; 1, 2];
for c = 1:9
    if c <= 3
        for r = 1:3
            S_macro(r,c)   = responses(c).stress_macro(r);
            S_macro(r+3,c) = sqrt(2) * responses(c).stress_macro(r+3);
            S_macro(r+6,c) = responses(c).E_macro(r);
        end
    elseif c <= 6
        for r = 1:3
            S_macro(r,c)   = sqrt(2) * responses(c).stress_macro(r);
            S_macro(r+3,c) = responses(c).stress_macro(r+3);
            S_macro(r+6,c) = responses(c).E_macro(r) / sqrt(2);
        end
    else
        for r = 1:3
            S_macro(r,c) = responses(c).stress_macro(r);
            S_macro(r+3,c) = sqrt(2) * responses(c).stress_macro(r+3);
            S_macro(r+6,c) = responses(c).E_macro(r);
        end
    end
end

% DERIVE MACROSCOPIC LINEAR CONSTITUTIVE MATRIX (ENTHALPY FORM)
S_mm = S_macro(1:6,1:6);
S_me = S_macro(1:6,7:9);
S_em = S_macro(7:9,1:6);
S_ee = S_macro(7:9,7:9);
C_mm = S_mm - (S_me / S_ee) * S_em;
C_me = S_me / S_ee;
C_em = - S_ee \ S_em;
C_ee = inv(S_ee);
C_macro = [C_mm, C_me; C_em, C_ee];

% fprintf('C_macro = \n'); fprintf([repmat('\t %8.5f ', 1, 9), '\n'], C_macro')
% fprintf('S_macro = \n'); fprintf([repmat('\t %8.5f ', 1, 9), '\n'], S_macro')
% file_name = sprintf('%06d_64x64x64.mat', responses(1).id);
file_name = sprintf('%s_64x64x64.mat', responses(1).id{1});
time = [responses(1).time, responses(2).time, responses(3).time, responses(4).time, responses(5).time, responses(6).time, responses(7).time, responses(8).time, responses(9).time];
% cd
% stress = responses(1).stress;
% E = responses(1).E;
cd('Results_test');
% save(fullfile('.', 'Results', file_name), 'C_macro', 'S_macro', 'time');
save(file_name, 'C_macro', 'S_macro', 'time', '-v7.3');
cd(current_path);
% save(fullfile('.', file_name), 'C_macro', 'S_macro', 'time', 'stress', 'E');
fprintf(log_file_id, '%s- Saved results \n', responses(1).id{1});

end
