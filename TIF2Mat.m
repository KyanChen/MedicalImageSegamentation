ROOTDIR = '/Users/keyanchen/Files/Ò½Ñ§Í¼Ïñ·Ö¸î/Segmentation_Result';

file_struct = dir(fullfile(ROOTDIR, '*.tif'));
Segmentation_Result = zeros(512, 512, length(file_struct), 'uint8');

for i=1:length(file_struct)
    img = imread([ROOTDIR, '/', file_struct(i).name]);
    Segmentation_Result(:,:,i) = img;
end
save('Segmentation_Result.mat');
imshow(Segmentation_Result(:,:,160));