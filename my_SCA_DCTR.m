function my_SCA_DCTR(image_path,feature_path,map_path, QF)
files=dir([image_path '\*.jpg']);
file_num=length(files)
names = cell(file_num,1);
beta= cell(file_num,1);
F = zeros(file_num,8000, 'single');

for w =1:file_num
    MAPfilename = [map_path '\' files(w).name(1:end-4) '.mat'];   %% jpegfilename = 'C:\Users\VAZ\Desktop\3.jpg';
    load(MAPfilename);   %%% 注意  这里 load的数据名称为Ori_beta；
    beta{w} = double(Ori_beta);
end

for w =1:file_num
    tic
    jpegfilename = [image_path '\' files(w).name];   %% jpegfilename = 'C:\Users\VAZ\Desktop\3.jpg';
    f = SCA_DCTR(jpegfilename, QF, beta{w});
    F(w,:) = f(:);   %将f以行向量的形式存储， w代表样本的个数。即每个样本的特征通过F对应一行的元素表示 %写成：的形式，就可以用parfor
    names{w} = files(w).name;
    toc
end

save(feature_path,'F','names','-v7.3');
disp('end')