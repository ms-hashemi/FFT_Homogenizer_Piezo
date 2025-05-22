function [phase_IND, varargout] = read_microimage_resize(path_to_images, prefix_image_files, padding_image_files, extension_image_files, boundaries, varargin)

%% Parsing the inputs
boundaries = int32(boundaries);
extracted_file_names = [];
n = -1;
[filepath, name, ext] = fileparts(path_to_images);
if strcmp(ext, '') % Means that if "path_to_images" is actually a folder containing the desired images
    % First check if the first slice can be found
    if isempty(dir(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files))))
        warning('MyComponent:notFoundFirstImage', ...
                [fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)) ' does not exist.']);
        phase_IND = -1; % Meaning that the function has encountered an error reading the images.
        if nargout > 1
            varargout{1} = n;
        end
        return
    % Now check if the last slice can be found
    elseif isempty(dir(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 2), extension_image_files))))
        warning('MyComponent:notFoundLastImage', ...
                [fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)) 'does not exist.']);
        phase_IND = -1; % Meaning that the function has encountered an error reading the images.
        if nargout > 1
            varargout{1} = n;
        end
        return
    % Now check if the requested boundaries are possible or within the size
    % of the images
    else
        image = imread(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)));
        image_size = size(image);
        if ((boundaries(1, 1) < 1) || (boundaries(1, 1) > image_size(1)) || ...
                (boundaries(1, 1) < 1) || (boundaries(1, 1) > image_size(1)) )
            warning('MyComponent:incorrectImageBoundaries', ...
                    'Image boundaries are out of bounds according to the first slice image.');
            phase_IND = -1; % Meaning that the function has encountered an error reading the images.
            if nargout > 1
                varargout{1} = n;
            end
            return
        end
    end
else % Meaning that "path_to_images" was not a folder, so try to extract the compressed file containing the images
    try
        extracted_file_names = unzip(path_to_images);
    catch
        try
            extracted_file_names = untar(path_to_images);
        catch ME
            warning('MyComponent:unableToExtract', ...
                    '%s', ME.message);
        end
    end
    if ~isempty(dir(fullfile(name)))
        path_to_images = fullfile(name);
    end
    % First check if the first slice can be found
    if isempty(dir(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files))))
        fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files))
        warning('MyComponent:notFoundFirstImage', ...
                [fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)) 'does not exist.']);
        phase_IND = -1; % Meaning that the function has encountered an error reading the images.
        if nargout > 1
            varargout{1} = n;
        end
        return
    % Now check if the last slice can be found
    elseif isempty(dir(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 2), extension_image_files))))
        warning('MyComponent:notFoundLastImage', ...
                [fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)) 'does not exist.']);
        phase_IND = -1; % Meaning that the function has encountered an error reading the images.
        if nargout > 1
            varargout{1} = n;
        end
        return
    % Now check if the requested boundaries are possible or within the size
    % of the images
    else
        image = imread(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, boundaries(3, 1), extension_image_files)));
        image_size = size(image);
        if ((boundaries(1, 1) < 1) || (boundaries(1, 1) > image_size(1)) || ...
                (boundaries(1, 1) < 1) || (boundaries(1, 1) > image_size(1)) )
            warning('MyComponent:incorrectImageBoundaries', ...
                    'Image boundaries are out of bounds according to the first slice image.');
            phase_IND = -1; % Meaning that the function has encountered an error reading the images.
            if nargout > 1
                varargout{1} = n;
            end
            return
        end
    end    
end

%% Clustering


%% Extracting the resultant matrix
% matrix phase is indicated by integer 1, IND stands for indices
phase_IND_original = zeros(boundaries(1, 2)-boundaries(1, 1)+1, boundaries(2, 2)-boundaries(2, 1)+1, boundaries(3, 2)-boundaries(3, 1)+1);
for k = boundaries(3, 1):1:boundaries(3, 2)
    image = imread(fullfile(path_to_images, sprintf('%s%0*d.%s', prefix_image_files, padding_image_files, k, extension_image_files)));
    for i = boundaries(1, 1):1:boundaries(1, 2)
        for j = boundaries(2, 1):1:boundaries(2, 2)
            if image(i, j) == 1
                phase_IND_original(i-boundaries(1, 1)+1, j-boundaries(2, 1)+1, k-boundaries(3, 1)+1) = 1;
            end
        end
    end
end
phase_IND_original = logical(phase_IND_original);

% Clean up the newly created files/folders if a compressed file has been
% extracted
% current_directory = cd;
if ~isempty(dir(path_to_images))
    rmdir(path_to_images, 's');
%     rmdir(fullfile(path_to_images, '..'), 's');
else
%     rmdir(path_to_images, 's');
end

%% Resize the image matrix

if nargin > 5
    n = varargin{1};
    phase_IND = imresize3(phase_IND_original, n, 'nearest');
else
    phase_IND = phase_IND_original;
    n = size(phase_IND);
end
phase_IND = int16(phase_IND);
phase_IND(phase_IND == 1) = 2;  % this line must preceed the below line as phase_IND contains 0 and 1 originally.
phase_IND(phase_IND == 0) = 1;

if nargout > 1
    varargout{1} = n;
end