% % Convert bmp to png
% parallel_size = 10;
% delete(gcp('nocreate'))
% p = parpool(parallel_size);
% mkdir('Vis');
% current_path = cd;
% 
% parfor bin = 1:number_bins
% %     current_path = cd;
%     for i_bin = 1:number_per_bin
%         boundaries = [16, 165; 16, 165; 16, 165];
%         if microstructure_folder_names(bin, i_bin) > 0
%             name = sprintf('%06d', microstructure_folder_names(bin, i_bin));
%         else
%             name = sprintf('_%06d', -microstructure_folder_names(bin, i_bin));
%         end
%         if i_bin >= 2 && microstructure_folder_names(bin, i_bin) == microstructure_folder_names(bin, i_bin-1)
%             continue
%         end
%         if ~isempty(dir(fullfile('Results', [name '.mat'])))
%             continue
%         end
%         unzip(fullfile('Results', name), 'Vis');
%         path_to_images = fullfile('Vis', name);
% %         path_to_images = fullfile('..', 'deg', 'Results', name);
%         prefix_image_files = 'slice_';
%         padding_image_files = 3;
%         extension_image_files = 'bmp';
%         phase_IND_original = zeros(boundaries(1, 2)-boundaries(1, 1)+1, boundaries(2, 2)-boundaries(2, 1)+1, boundaries(3, 2)-boundaries(3, 1)+1);
%         
%         if isempty(dir(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 2), extension_image_files))))
%             boundaries = [1, 150; 1, 150; 1, 150];
%             extension_image_files = 'png';
%             for k = boundaries(3, 1):1:boundaries(3, 2)
%                 image = imread(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, k, extension_image_files)));
%                 for i = boundaries(1, 1):1:boundaries(1, 2)
%                     for j = boundaries(2, 1):1:boundaries(2, 2)
%                         if image(i, j) == 1
%                             phase_IND_original(i-boundaries(1, 1)+1, j-boundaries(2, 1)+1, k-boundaries(3, 1)+1) = 1;
%                         end
%                     end
%                 end
%             end
%             phase_IND_original = logical(phase_IND_original);
%             save_results(phase_IND_original, name, 'Results', current_path);
%             rmdir(path_to_images, 's');
%             continue
%         end
% 
%         for k = boundaries(3, 1):1:boundaries(3, 2)
%             image = imread(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, k, extension_image_files)));
%             for i = boundaries(1, 1):1:boundaries(1, 2)
%                 for j = boundaries(2, 1):1:boundaries(2, 2)
%                     if image(i, j) == 1
%                         phase_IND_original(i-boundaries(1, 1)+1, j-boundaries(2, 1)+1, k-boundaries(3, 1)+1) = 1;
%                     end
%                 end
%             end
%         end
%         phase_IND_original = logical(phase_IND_original);
%         
%         VF = sum(phase_IND_original, 'all')/(150^3);
%         if VF < 0.2 || VF > 0.9
%             fprintf('%s: %.2f', name, VF);
%         end
%         
%         path_to_results = fullfile('Results', name);
%         mkdir(path_to_results);
%         map = [0, 0, 0; 1, 1, 1];
%         extension_image_files = 'png';
%         for k = 1:1:boundaries(3, 2)-boundaries(3, 1)+1
%             imwrite(phase_IND_original(:, :, k), map, fullfile(path_to_results, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, k, extension_image_files)));
%         end
%         zip([path_to_results '.zip'], name, 'Results');
%         save_results(phase_IND_original, name, 'Results', current_path);
%         rmdir(path_to_results, 's');
%         rmdir(path_to_images, 's');
%     end
% end
% 
% function save_results(results, file_name, save_path, current_path)
%     cd(save_path);
%     phase_IND_original = results;
%     save(file_name, 'phase_IND_original', '-v7.3');
%     cd(current_path);
% end

% % Stress contour
% load('000011_64x64x64.mat');
% imshow(logical(phase_IND(:,:,32)-1));
% stress = reshape(stress, [[64, 64, 64], 6]);
% [X, Y] = meshgrid(1:64, 1:64);
% figure;
% axis equal
% contourf(X,Y,stress(:, :, 32, 1), 10);

% % Calculating propery ranges
% outputs_path = fullfile('Results');
% folder_info = dir(fullfile(outputs_path, '*_64x64x64.mat'));
% for i = 1:1:length(folder_info)
%     variables = load(fullfile(outputs_path, folder_info(i).name), 'C_macro');
%     C11(i) = (variables.C_macro(1, 1) + variables.C_macro(2, 2))/2;
%     C12(i) = (variables.C_macro(1, 2) + variables.C_macro(2, 1))/2;
%     C13(i) = (variables.C_macro(1, 3) + variables.C_macro(3, 1) + variables.C_macro(2, 3) + variables.C_macro(3, 2))/4;
%     C33(i) = variables.C_macro(3, 3);
%     C44(i) = (variables.C_macro(4, 4) + variables.C_macro(5, 5))/2;
%     C66(i) = variables.C_macro(6, 6);
%     e31(i) = (variables.C_macro(9, 1) + variables.C_macro(9, 2))/2;
%     e33(i) = variables.C_macro(9, 3);
%     e15(i) = (variables.C_macro(8, 4) + variables.C_macro(7, 5))/2;
%     gamma11(i) = (variables.C_macro(7, 7) + variables.C_macro(8, 8))/2;
%     gamma33(i) = variables.C_macro(9, 9);
% %     clear('variables');
% ends
% [X,Y] = meshgrid(min(C33):1:max(C33), min(C11):1:max(C11));
% steps_x = length(min(C33):1:max(C33)) - 1;
% steps_y = length(min(C11):1:max(C11)) - 1;
% % Plotting the Z-values of the function whose level sets have
% % to be determined
% Z = zeros(steps_y, steps_x);
% for i = 1:1:length(folder_info)
%     for j = 1:1:steps_x
%         for k = 1:1:steps_y
%             if C33(i) >= X(1, j) && C33(i) < X(1, j+1) && C11(i) >= Y(k, 1) && C11(i) < Y(k+1, 1)
%                 Z(k, j) = Z(k, j) + 1;
%             end
%         end
%     end
% end
% pcolor(X(1:end-1,1:end-1), Y(1:end-1,1:end-1), Z)
% colorbar

% n = [150, 150, 150];
% load('test.mat');
% strain = reshape(strain, [n, 6]);
% D = reshape(D, [n, 3]);
% stress = reshape(stress, [n, 6]);
% E = reshape(E, [n, 3]);
% % save('test2.mat', 'stress', 'strain', 'D', 'E');
% 
% x = 0:1/n(1):1;
% y = 0:1/n(2):1;
% z = 0:1/n(3):1;
% [X, Y, Z] = meshgrid(x, y, z);
% % vtkwrite('test.vtk', 'structured_grid', X, Y, Z, 'scalars', 'S11', stress(:, :, :, 1), 'S22', stress(:, :, :, 2), 'S33', stress(:, :, :, 3), 'S23', stress(:, :, :, 4), 'S13', stress(:, :, :, 5), 'S12', stress(:, :, :, 6)...
% %     , 'Strain11', strain(:, :, :, 1), 'Strain22', strain(:, :, :, 2), 'Strain33', strain(:, :, :, 3), 'Strain23', strain(:, :, :, 4), 'Strain13', strain(:, :, :, 5), 'Strain12', strain(:, :, :, 6)...
% %     , 'vectors', 'E', E(:, :, :, 1), E(:, :, :, 2), E(:, :, :, 3), 'D', D(:, :, :, 1), D(:, :, :, 2), D(:, :, :, 3));
% 
% DT = delaunayTriangulation(X(:), Y(:), Z(:));
% % vtkwrite('test2.vtk', 'polydata', 'tetrahedron', X, Y, Z, DT.ConnectivityList, 'scalars', 'S11', stress(:, :, :, 1), 'scalars', 'S22', stress(:, :, :, 2), 'scalars', 'S33', stress(:, :, :, 3), 'scalars', 'S23', stress(:, :, :, 4), 'scalars', 'S13', stress(:, :, :, 5), 'scalars', 'S12', stress(:, :, :, 6)...
% %     , 'scalars', 'Strain11', strain(:, :, :, 1), 'scalars', 'Strain22', strain(:, :, :, 2), 'scalars', 'Strain33', strain(:, :, :, 3), 'scalars', 'Strain23', strain(:, :, :, 4), 'scalars', 'Strain13', strain(:, :, :, 5), 'scalars', 'Strain12', strain(:, :, :, 6)...
% %     , 'vectors', 'E', E(:, :, :, 1), E(:, :, :, 2), E(:, :, :, 3), 'vectors', 'D', D(:, :, :, 1), D(:, :, :, 2), D(:, :, :, 3));
% 
% % vtkwrite('test3.vtk', 'polydata', 'tetrahedron', X, Y, Z, DT.ConnectivityList, 'structured_grid', 'scalars', 'S11', stress(:, :, :, 1));
% 
% % vtkwrite('test.vtk', 'structured_points', 'Stress11', stress(:, :, :, 1), 'spacing', sx, sy, sz);
% 
% % cmap = colormap('jet');
% % vol = volshow(255*(stress(:, :, :, 1)-min(reshape(stress(:, :, :, 1), [150^3, 1])))/(max(reshape(stress(:, :, :, 1), [150^3, 1]))-min(reshape(stress(:, :, :, 1), [150^3, 1]))), 'Colormap', colormap);
% 
% data = zeros(n(1)*n(2)*n(3), 5);
% for i = 1:1:n(1)
%     for j = 1:1:n(2)
%         for k = 1:1:n(3)
%             data((i-1)*(n(2)*n(3)) + (j-1)*n(3) + k, :) = [i, j, k, stress(i, j, k, 1), phase_IND(i, j, k)-1];
%         end
%     end
% end
% 
% 
% writematrix(data, 'data.csv');
% 
% data2 = zeros(n(1)*n(2)*n(3), 5);
% for i = 1:1:n(1)
%     for j = 1:1:n(2)
%         for k = 1:1:n(3)
%             data2((i-1)*(n(2)*n(3)) + (j-1)*n(3) + k, :) = [i, j, k, E(i, j, k, 1), phase_IND(i, j, k)-1];
%         end
%     end
% end
% 
% 
% writematrix(data2, 'data2.csv');

% Training results
path_to_log_files = '../../model/TVG_Design[64, 64, 64]/TVG_Conv-ViT-Gen2-B_16_vitpatch[1, 1, 1]_epo200_bs35_lr0.001_seed1234';
file_id = fopen(fullfile(path_to_log_files, 'log.txt'), 'r');
first_line = 1;    
number_of_trainings = 0;
while ~feof(file_id) && ~all(first_line == -1)
    before_first_line = fgets(file_id);
    position = ftell(file_id);
    first_line = fgets(file_id);
    fseek(file_id, position, 'bof');
    if ~all(first_line == -1)
        extracted_data = textscan([before_first_line first_line], '[%d:%d:%f] iteration %d: loss: %f, loss_kl: %f, loss_recon: %f, loss_pred: %f\n[%d:%d:%f] iteration %d: loss: %*f, loss_kl: %*f, loss_recon: %*f, loss_pred: %*f\n');
        if ~any(cellfun(@isempty, extracted_data))
            if number_of_trainings == 0
                number_of_trainings = number_of_trainings + 1;
                times(number_of_trainings) = 0;
                loss(number_of_trainings, extracted_data{4}) = extracted_data{5};
                loss_kl(number_of_trainings, extracted_data{4}) = extracted_data{6};
                loss_recon(number_of_trainings, extracted_data{4}) = extracted_data{7};
                loss_pred(number_of_trainings, extracted_data{4}) = extracted_data{8};
            elseif extracted_data{4} ~= extracted_data{12}-1
                loss(number_of_trainings, extracted_data{4}) = extracted_data{5};
                loss_kl(number_of_trainings, extracted_data{4}) = extracted_data{6};
                loss_recon(number_of_trainings, extracted_data{4}) = extracted_data{7};
                loss_pred(number_of_trainings, extracted_data{4}) = extracted_data{8};
                number_of_trainings = number_of_trainings + 1;
                times(number_of_trainings) = 0;
            else
                if extracted_data{9} < extracted_data{1}
                    times(number_of_trainings) = (extracted_data{9}*3600 + extracted_data{10}*60 + extracted_data{11}) + (86400 - extracted_data{1}*3600 + extracted_data{2}*60 + extracted_data{3}) + times(number_of_trainings);
                else
                    times(number_of_trainings) = (extracted_data{9}*3600 + extracted_data{10}*60 + extracted_data{11}) - (extracted_data{1}*3600 + extracted_data{2}*60 + extracted_data{3}) + times(number_of_trainings);
                end
                loss(number_of_trainings, extracted_data{4}) = extracted_data{5};
                loss_kl(number_of_trainings, extracted_data{4}) = extracted_data{6};
                loss_recon(number_of_trainings, extracted_data{4}) = extracted_data{7};
                loss_pred(number_of_trainings, extracted_data{4}) = extracted_data{8};
            end
        else
            extracted_data = textscan(before_first_line, '[%d:%d:%f] iteration %d: loss: %f, loss_kl: %f, loss_recon: %f, loss_pred: %f\n');
            if ~any(cellfun(@isempty, extracted_data))
                if number_of_trainings == 0
                    number_of_trainings = number_of_trainings + 1;
                    times(number_of_trainings) = 0;
                    loss(number_of_trainings, extracted_data{4}) = extracted_data{5};
                    loss_kl(number_of_trainings, extracted_data{4}) = extracted_data{6};
                    loss_recon(number_of_trainings, extracted_data{4}) = extracted_data{7};
                    loss_pred(number_of_trainings, extracted_data{4}) = extracted_data{8};
                    number_of_trainings = number_of_trainings + 1;
                else
                    loss(number_of_trainings, extracted_data{4}) = extracted_data{5};
                    loss_kl(number_of_trainings, extracted_data{4}) = extracted_data{6};
                    loss_recon(number_of_trainings, extracted_data{4}) = extracted_data{7};
                    loss_pred(number_of_trainings, extracted_data{4}) = extracted_data{8};
                end
            end
        end
    end
end
fclose(file_id);

max_epochs = 200;
start = 0.01;
stop = 1.0;
n_cycle = 4;
ratio = 0.5;
r = ones(1, max_epochs);
period = floor(max_epochs / n_cycle);
step = (stop - start) / (period * ratio); % step is in [0, 1]
% Sigmoid cyclical
% transform into [-6, 6] for plots: v*12.-6.
for c = 0:1:n_cycle-1
    v = start;
    i = 1;
    while v <= stop
        r(floor(i + c*period)) = 1.0/(1.0+ exp(- (v*12.-6.)));
        v = v + step;
        i = i + 1;
    end
end
% Linear cyclical
for c = 0:1:n_cycle-1
    v = start;
    i = 1;
    while v <= stop
        if i < ((1 - ratio)/2) * period - 1
            r(floor(i+c*period)) = start;
        else
            r(floor(i+c*period)) = v;
            v = v + step;
        end
        i = i + 1;
    end
end
figure;
plot(r, '-k', 'LineWidth', 2);
xlabel('Epoch');
ylabel('r');

for trial = 1:1:number_of_trainings
    figure;
    sz = size(loss);
    iterations = 1:1:sz(2);
    plot(iterations, loss(trial, :), iterations, 10*loss_kl(trial, :), iterations, -loss_recon(trial, :), iterations, 100*loss_pred(trial, :));
    xlabel('Iteration');
    ylabel('Loss');
    xlim([0 sz(2)])
    ylim([0 inf]);
    legend('Total', 'KL', 'Reconstruction', 'Prediction')
end

% % Regression plots
% path_to_log_files = '../../TransVNet/test_log/TVG_Design[64, 64, 64]';
% file_id = fopen(fullfile(path_to_log_files, 'TVG_Conv-ViT-Gen2-B_16_vitpatch[1, 1, 1]_epo200_bs35_lr0.001_seed1234.txt'), 'r');
% % file_id = fopen(fullfile(path_to_log_files, 'new 1.txt'), 'r');
% first_line = 1;
% counter = 0;
% name = {};
% while ~feof(file_id) && ~all(first_line == -1)
%     before_first_line = fgets(file_id);
%     position = ftell(file_id);
%     first_line = fgets(file_id);
%     fseek(file_id, position, 'bof');
%     counter_before = counter;
%     if ~all(first_line == -1)
%         extracted_data = textscan([before_first_line first_line], '[%*d:%*d:%*f] name %s %f surrogate_model_error %f generative_error %f reconstruction_loss %f C11 %f C11_surrogate %f C11_surrogate_error %f C11_generative %f C11_generative_error %f C12 %f C12_surrogate %f C12_surrogate_error %f C12_generative %f C12_generative_error %f C13 %f C13_surrogate %f C13_surrogate_error %f C13_generative %f C13_generative_error %f C33 %f C33_surrogate %f C33_surrogate_error %f C33_generative %f C33_generative_error %f C44 %f C44_surrogate %f C44_surrogate_error %f C44_generative %f C44_generative_error %f C66 %f C66_surrogate %f C66_surrogate_error %f C66_generative %f C66_generative_error %f e31 %f e31_surrogate %f e31_surrogate_error %f e31_generative %f e31_generative_error %f e33 %f e33_surrogate %f e33_surrogate_error %f e33_generative %f e33_generative_error %f e15 %f e15_surrogate %f e15_surrogate_error %f e15_generative %f e15_generative_error %f gamma11 %f gamma11_surrogate %f gamma11_surrogate_error %f gamma11_generative %f gamma11_generative_error %f gamma33 %f gamma33_surrogate %f gamma33_surrogate_error %f gamma33_generative %f gamma33_generative_error %f\n[%*d:%*d:%*f] name %s %f surrogate_model_error %f generative_error %f reconstruction_loss %f C11 %f C11_surrogate %f C11_surrogate_error %f C11_generative %f C11_generative_error %f C12 %f C12_surrogate %f C12_surrogate_error %f C12_generative %f C12_generative_error %f C13 %f C13_surrogate %f C13_surrogate_error %f C13_generative %f C13_generative_error %f C33 %f C33_surrogate %f C33_surrogate_error %f C33_generative %f C33_generative_error %f C44 %f C44_surrogate %f C44_surrogate_error %f C44_generative %f C44_generative_error %f C66 %f C66_surrogate %f C66_surrogate_error %f C66_generative %f C66_generative_error %f e31 %f e31_surrogate %f e31_surrogate_error %f e31_generative %f e31_generative_error %f e33 %f e33_surrogate %f e33_surrogate_error %f e33_generative %f e33_generative_error %f e15 %f e15_surrogate %f e15_surrogate_error %f e15_generative %f e15_generative_error %f gamma11 %f gamma11_surrogate %f gamma11_surrogate_error %f gamma11_generative %f gamma11_generative_error %f gamma33 %f gamma33_surrogate %f gamma33_surrogate_error %f gamma33_generative %f gamma33_generative_error %f\n');
%         if ~any(cellfun(@isempty, extracted_data))
%             counter = counter + 1;
%             data.name(counter) = extracted_data{1};
%             data.VF(counter) = extracted_data{2};
%             data.surrogate_model_error(counter) = extracted_data{3};
%             data.generative_error(counter) = extracted_data{4};
%             data.reconstruction_loss(counter) = extracted_data{5};
%             data.C11(counter) = extracted_data{6};
%             data.C11_surrogate(counter) = extracted_data{7};
%             data.C11_surrogate_error(counter) = extracted_data{8};
%             data.C11_generative(counter) = extracted_data{9};
%             data.C11_generative_error(counter) = extracted_data{10};
%             data.C12(counter) = extracted_data{11};
%             data.C12_surrogate(counter) = extracted_data{12};
%             data.C12_surrogate_error(counter) = extracted_data{13};
%             data.C12_generative(counter) = extracted_data{14};
%             data.C12_generative_error(counter) = extracted_data{15};
%             data.C13(counter) = extracted_data{16};
%             data.C13_surrogate(counter) = extracted_data{17};
%             data.C13_surrogate_error(counter) = extracted_data{18};
%             data.C13_generative(counter) = extracted_data{19};
%             data.C13_generative_error(counter) = extracted_data{20};
%             data.C33(counter) = extracted_data{21};
%             data.C33_surrogate(counter) = extracted_data{22};
%             data.C33_surrogate_error(counter) = extracted_data{23};
%             data.C33_generative(counter) = extracted_data{24};
%             data.C33_generative_error(counter) = extracted_data{25};
%             data.C44(counter) = extracted_data{26};
%             data.C44_surrogate(counter) = extracted_data{27};
%             data.C44_surrogate_error(counter) = extracted_data{28};
%             data.C44_generative(counter) = extracted_data{29};
%             data.C44_generative_error(counter) = extracted_data{30};
%             data.C66(counter) = extracted_data{31};
%             data.C66_surrogate(counter) = extracted_data{32};
%             data.C66_surrogate_error(counter) = extracted_data{33};
%             data.C66_generative(counter) = extracted_data{34};
%             data.C66_generative_error(counter) = extracted_data{35};
%             data.e31(counter) = extracted_data{36};
%             data.e31_surrogate(counter) = extracted_data{37};
%             data.e31_surrogate_error(counter) = extracted_data{38};
%             data.e31_generative(counter) = extracted_data{39};
%             data.e31_generative_error(counter) = extracted_data{40};
%             data.e33(counter) = extracted_data{41};
%             data.e33_surrogate(counter) = extracted_data{42};
%             data.e33_surrogate_error(counter) = extracted_data{43};
%             data.e33_generative(counter) = extracted_data{44};
%             data.e33_generative_error(counter) = extracted_data{45};
%             data.e15(counter) = extracted_data{46};
%             data.e15_surrogate(counter) = extracted_data{47};
%             data.e15_surrogate_error(counter) = extracted_data{48};
%             data.e15_generative(counter) = extracted_data{49};
%             data.e15_generative_error(counter) = extracted_data{50};
%             data.gamma11(counter) = extracted_data{51};
%             data.gamma11_surrogate(counter) = extracted_data{52};
%             data.gamma11_surrogate_error(counter) = extracted_data{53};
%             data.gamma11_generative(counter) = extracted_data{54};
%             data.gamma11_generative_error(counter) = extracted_data{55};
%             data.gamma33(counter) = extracted_data{56};
%             data.gamma33_surrogate(counter) = extracted_data{57};
%             data.gamma33_surrogate_error(counter) = extracted_data{58};
%             data.gamma33_generative(counter) = extracted_data{59};
%             data.gamma33_generative_error(counter) = extracted_data{60};
%         else
%             extracted_data = textscan(before_first_line, '[%*d:%*d:%*f] name %s %f surrogate_model_error %f generative_error %f reconstruction_loss %f C11 %f C11_surrogate %f C11_surrogate_error %f C11_generative %f C11_generative_error %f C12 %f C12_surrogate %f C12_surrogate_error %f C12_generative %f C12_generative_error %f C13 %f C13_surrogate %f C13_surrogate_error %f C13_generative %f C13_generative_error %f C33 %f C33_surrogate %f C33_surrogate_error %f C33_generative %f C33_generative_error %f C44 %f C44_surrogate %f C44_surrogate_error %f C44_generative %f C44_generative_error %f C66 %f C66_surrogate %f C66_surrogate_error %f C66_generative %f C66_generative_error %f e31 %f e31_surrogate %f e31_surrogate_error %f e31_generative %f e31_generative_error %f e33 %f e33_surrogate %f e33_surrogate_error %f e33_generative %f e33_generative_error %f e15 %f e15_surrogate %f e15_surrogate_error %f e15_generative %f e15_generative_error %f gamma11 %f gamma11_surrogate %f gamma11_surrogate_error %f gamma11_generative %f gamma11_generative_error %f gamma33 %f gamma33_surrogate %f gamma33_surrogate_error %f gamma33_generative %f gamma33_generative_error %f\n');
%             if ~any(cellfun(@isempty, extracted_data))
%                 counter = counter + 1;
%                 data.name(counter) = extracted_data{1};
%                 data.VF(counter) = extracted_data{2};
%                 data.surrogate_model_error(counter) = extracted_data{3};
%                 data.generative_error(counter) = extracted_data{4};
%                 data.reconstruction_loss(counter) = extracted_data{5};
%                 data.C11(counter) = extracted_data{6};
%                 data.C11_surrogate(counter) = extracted_data{7};
%                 data.C11_surrogate_error(counter) = extracted_data{8};
%                 data.C11_generative(counter) = extracted_data{9};
%                 data.C11_generative_error(counter) = extracted_data{10};
%                 data.C12(counter) = extracted_data{11};
%                 data.C12_surrogate(counter) = extracted_data{12};
%                 data.C12_surrogate_error(counter) = extracted_data{13};
%                 data.C12_generative(counter) = extracted_data{14};
%                 data.C12_generative_error(counter) = extracted_data{15};
%                 data.C13(counter) = extracted_data{16};
%                 data.C13_surrogate(counter) = extracted_data{17};
%                 data.C13_surrogate_error(counter) = extracted_data{18};
%                 data.C13_generative(counter) = extracted_data{19};
%                 data.C13_generative_error(counter) = extracted_data{20};
%                 data.C33(counter) = extracted_data{21};
%                 data.C33_surrogate(counter) = extracted_data{22};
%                 data.C33_surrogate_error(counter) = extracted_data{23};
%                 data.C33_generative(counter) = extracted_data{24};
%                 data.C33_generative_error(counter) = extracted_data{25};
%                 data.C44(counter) = extracted_data{26};
%                 data.C44_surrogate(counter) = extracted_data{27};
%                 data.C44_surrogate_error(counter) = extracted_data{28};
%                 data.C44_generative(counter) = extracted_data{29};
%                 data.C44_generative_error(counter) = extracted_data{30};
%                 data.C66(counter) = extracted_data{31};
%                 data.C66_surrogate(counter) = extracted_data{32};
%                 data.C66_surrogate_error(counter) = extracted_data{33};
%                 data.C66_generative(counter) = extracted_data{34};
%                 data.C66_generative_error(counter) = extracted_data{35};
%                 data.e31(counter) = extracted_data{36};
%                 data.e31_surrogate(counter) = extracted_data{37};
%                 data.e31_surrogate_error(counter) = extracted_data{38};
%                 data.e31_generative(counter) = extracted_data{39};
%                 data.e31_generative_error(counter) = extracted_data{40};
%                 data.e33(counter) = extracted_data{41};
%                 data.e33_surrogate(counter) = extracted_data{42};
%                 data.e33_surrogate_error(counter) = extracted_data{43};
%                 data.e33_generative(counter) = extracted_data{44};
%                 data.e33_generative_error(counter) = extracted_data{45};
%                 data.e15(counter) = extracted_data{46};
%                 data.e15_surrogate(counter) = extracted_data{47};
%                 data.e15_surrogate_error(counter) = extracted_data{48};
%                 data.e15_generative(counter) = extracted_data{49};
%                 data.e15_generative_error(counter) = extracted_data{50};
%                 data.gamma11(counter) = extracted_data{51};
%                 data.gamma11_surrogate(counter) = extracted_data{52};
%                 data.gamma11_surrogate_error(counter) = extracted_data{53};
%                 data.gamma11_generative(counter) = extracted_data{54};
%                 data.gamma11_generative_error(counter) = extracted_data{55};
%                 data.gamma33(counter) = extracted_data{56};
%                 data.gamma33_surrogate(counter) = extracted_data{57};
%                 data.gamma33_surrogate_error(counter) = extracted_data{58};
%                 data.gamma33_generative(counter) = extracted_data{59};
%                 data.gamma33_generative_error(counter) = extracted_data{60};
%             end
%         end
%         if counter_before ~= counter - 1
%             counter = 0;
%             data.name = {};
%             data.VF = [];
%             data.surrogate_model_error = [];
%             data.generative_error = [];
%             data.reconstruction_loss = [];
%             data.C11 = [];
%             data.C11_surrogate = [];
%             data.C11_surrogate_error = [];
%             data.C11_generative = [];
%             data.C11_generative_error = [];
%             data.C12 = [];
%             data.C12_surrogate = [];
%             data.C12_surrogate_error = [];
%             data.C12_generative = [];
%             data.C12_generative_error = [];
%             data.C13 = [];
%             data.C13_surrogate = [];
%             data.C13_surrogate_error = [];
%             data.C13_generative = [];
%             data.C13_generative_error = [];
%             data.C33 = [];
%             data.C33_surrogate = [];
%             data.C33_surrogate_error = [];
%             data.C33_generative = [];
%             data.C33_generative_error = [];
%             data.C44 = [];
%             data.C44_surrogate = [];
%             data.C44_surrogate_error = [];
%             data.C44_generative = [];
%             data.C44_generative_error = [];
%             data.C66 = [];
%             data.C66_surrogate = [];
%             data.C66_surrogate_error = [];
%             data.C66_generative = [];
%             data.C66_generative_error = [];
%             data.e31 = [];
%             data.e31_surrogate = [];
%             data.e31_surrogate_error = [];
%             data.e31_generative = [];
%             data.e31_generative_error = [];
%             data.e33 = [];
%             data.e33_surrogate = [];
%             data.e33_surrogate_error = [];
%             data.e33_generative = [];
%             data.e33_generative_error = [];
%             data.e15 = [];
%             data.e15_surrogate = [];
%             data.e15_surrogate_error = [];
%             data.e15_generative = [];
%             data.e15_generative_error = [];
%             data.gamma11 = [];
%             data.gamma11_surrogate = [];
%             data.gamma11_surrogate_error = [];
%             data.gamma11_generative = [];
%             data.gamma11_generative_error = [];
%             data.gamma33 = [];
%             data.gamma33_surrogate = [];
%             data.gamma33_surrogate_error = [];
%             data.gamma33_generative = [];
%             data.gamma33_generative_error = [];
%         end
%     end
% end
% fclose(file_id);
% 
% % Show the performance of the generative network on the test dataset
% % desired_data_fields = fieldnames(data);
% properties = {'C33', 'e33'};
% model = {'_surrogate', '_generative'};
% 
% for i = 1:1:numel(properties)
%     for j = 1:1:numel(model)
%         x = (data.(properties{i}))';
%         y = (data.([properties{i} model{j}]))';
%         % Scatter plot
%         figure;
%         scatter(x, y, 40, 100*(data.VF)', 'filled', 's', 'LineWidth', 3);
%         ax = gca;
%         ax.FontSize = 14;
%         ax.FontName = 'Times New Roman';
%         xlabel(['FFT-Homogenized ' properties{i}]);
%         if strcmp(model{j}, '_surrogate')
%             ylabel(['Surrogate-Predicted ' properties{i}]);
%         else
%             ylabel(['Generative-Estimated ' properties{i}]);
%         end
%         axis([0 ceil(max([x' y'])) 0 ceil(max([x' y']))]);
%         axis tight
%         axis equal
%         colormap('gray');
%         c = colorbar('eastoutside');
%         c.Label.String = 'Volume Fraction in the test dataset (%)';
%         % c.Label.FontSize = 12;
% 
%         % Fit linear regression line with OLS.
%         % x = linspace(0, ceil(max(C33)), 1001);
%         % b = [ones(size(x,1), 1) x]\y;
%         b = [ones(size(x,1), 1) x]\x;
%         % Use estimated slope and intercept to create regression line.
%         RegressionLine = [ones(size(x,1), 1) x]*b;
%         % RegressionLine = x;
%         % Plot it in the scatter plot and show equation.
%         hold on
%         if abs(b(1)) < 1000*eps
%             p = plot(x, RegressionLine, 'k', 'linewidth', 3, 'displayname', sprintf('Regression line (y = x)'));
%         %     p = plot(x, RegressionLine, 'displayname', sprintf('Regression line (y = %.3f*x)', b(2)));
%         else
%             p = plot(x, RegressionLine, 'displayname', sprintf('Regression line (y = %.3f*x + %.3f)', b(2), b(1)));
%         end
%         leg = legend('Observed values', p.DisplayName);%, '95% C.I'
%         % set(leg, 'location', 'best')
%         set(leg, 'location', 'southeast')
%         % RMSE between regression line and y
%         RMSE = sqrt(mean((y - RegressionLine).^2));
%         % R2 between regression line and y
%         SS_X = sum((RegressionLine - mean(RegressionLine)).^2);
%         SS_Y = sum((y - mean(y)).^2);
%         SS_XY = sum((RegressionLine - mean(RegressionLine)).*(y - mean(y)));
%         R_squared = SS_XY / sqrt(SS_X * SS_Y);
%         fprintf('RMSE: %0.2f | R2: %0.2f\n', RMSE, R_squared)
%         bias = sum(y - x)/length(y);
%         str = ['N = ', sprintf('%d', length(x)),...
%             ', Bias = ', sprintf('%.3f', bias),...    
%             ', R^2 = ', sprintf('%.2f', R_squared),...
%             ', RMSE = ', sprintf('%.2f', RMSE),
%             %'y = ',sprintf('%.2f',table2array(mdl.Coefficients(2,1))),'x + ',sprintf('%.2f',table2array(mdl.Coefficients(1,1)))
%         ]
%         annotation('textbox', [.2 0.9 0 0], 'string', str, 'FitBoxToText', 'on', 'EdgeColor', 'black', 'FontSize', 12, 'FontName', 'Times New Roman')
% 
%         figure;
%         X = [ones(size(x)) x];
%         [b, bint] = regress(y, X)
%         xval = min(x) - 0.05:0.01:max(x) + 0.05;
%         yhat = b(1) + b(2) * xval;
%         ylow = bint(1,1) + bint(2,1) * xval;
%         yupp = bint(1,2) + bint(2,2) * xval;
%         plot(x, y, 'ks', 'LineWidth', 3, 'MarkerSize', 2);
%         ax = gca;
%         ax.FontSize = 14;
%         ax.FontName = 'Times New Roman';
%         % p5.Color(4) = 0.5;
%         hold on;
%         p6 = plot(xval, yhat, 'k', 'linewidth', 3);
%         p6.Color(4) = 0.5;
%         % axis([0.04 0.3 0.03 .35])
%         % fontSize = 12;
%         hold on
%         patch([xval fliplr(xval)], [ylow fliplr(yupp)], 'k', 'EdgeColor', 'white')
%         alpha(0.3)
%         leg = legend('Observed values','Regression line','95% C.I');
%         set(leg, 'location', 'best')
%         xlabel(['Target ' properties{i}]);
%         if strcmp(model{j}, '_surrogate')
%             ylabel(['Surrogate-Predicted ' properties{i}]);
%         else
%             ylabel(['Generative-Estimated ' properties{i}]);
%         end
%         set(gcf, 'color', 'white')
%         bias = sum(y - x)/length(y);
%         tbl = table(y, x)
%         mdl = fitlm(tbl,'linear')
%         str = [    'N = ',sprintf('%d',mdl.NumObservations),...
%         ', Bias = ',sprintf('%.3f',bias),...    
%         ', R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary),...
%         %'y = ',sprintf('%.2f',table2array(mdl.Coefficients(2,1))),'x + ',sprintf('%.2f',table2array(mdl.Coefficients(1,1)))
%         ]
%         annotation('textbox', [.2 0.9 0 0], 'string', str, 'FitBoxToText', 'on', 'EdgeColor', 'black')
%     end
% end
