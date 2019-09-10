% test
clc;
clear;
% DICOM文件所在路径
ROOTPATH = 'E:\DataSet\医学图像分割\WWZ\IWOB2LVA\PV1YAB2C';

% 图片名字的正则表达式，需要修改
% 匹配 I9900000 类型的
Pattern = 'I\d{7}';

% 匹配image_29540375455130352.dcm类型
% Pattern = 'image_\d+.dcm';

imgPath = [ROOTPATH, '\TIF'];
if ~exist(imgPath, 'dir')
    mkdir(imgPath);
end
numImg = 1;
filelist = dir(ROOTPATH);
for i=1:length(filelist)
    str = filelist(i).name;
    matchStr = regexp(str, Pattern, 'match');
    % 如果匹配不成功则匹配下一个
    if isempty(matchStr)
        continue;
    end
    % 如果匹配成功，而且字符个数匹配上，则进行图片处理
    if length(matchStr{1}) == length(str)
        img=dicomread([ROOTPATH, '\', matchStr{1}]);
        img(img<0) = 0;
        img(img>2048) = 2047;
        img = uint16(img*2^3);
        fileNum = sprintf('%03d', numImg);
        filePath = [imgPath, '\result', fileNum, '.tif'];
        imwrite(img, filePath, 'Compression','none');
        numImg = numImg + 1;
    end
end