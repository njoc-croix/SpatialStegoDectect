clc;
clear;

% Set the paths for cover images and stego images folders
Input_cover_images = 'C:\Users\JEANDELACROIXNTIVUGU\Desktop\SIMPA for Cyber\Covers\';
Output_stego_images = 'C:\Users\JEANDELACROIXNTIVUGU\Desktop\SIMPA for Cyber\Stego\';

% Get a list of cover image files in the cover folder
cover_images = dir(fullfile(Input_cover_images, '*.pgm'));

% Set payload
payload_size = 0.4;

% Set params
params.p = -1;  % Holder norm parameter

fprintf('Random bitstreams Embedding using MATLAB\n');

for i = 1:numel(cover_images)
    coverPath = fullfile(Input_cover_images, cover_images(i).name);
    cover = imread(coverPath);
    
    fprintf('Data embedding in an inquiry image %d\n', i);
    
    % Run embedding simulation
    [stego, distortion] = WOW_main_function(cover, payload_size, params);
    
    % Save the stego image in the stego images folder with a modified name
    [~, filename, ext] = fileparts(cover_images(i).name);
    stego_image_name = [filename, '_stego', ext];
    stegoPath = fullfile(Output_stego_images, stego_image_name);
    imwrite(stego, stegoPath);
    
    fprintf('Stego image %d saved\n', i);
    
    % Display cover and stego images
    figure;
    subplot(1, 2, 1);
    imshow(cover, []);
    title('Cover');
    
    subplot(1, 2, 2);
    imshow(stego, []);
    title('Stego image');
    
    % Display changes resulted from embedding
    figure;
    imshow((double(stego) - double(cover) + 1)/2);
    title('Changes Resulted from Embedding');
    
    fprintf('Image %d embedded, change rate: %.4f, distortion per pixel: %.6f\n', i, sum(cover(:) ~= stego(:))/numel(cover), distortion/numel(cover));
end

fprintf('All images embedded using MATLAB code\n');
