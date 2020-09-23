function my_QF_GFR(image_path,feature_path, QF, num_size)
files=dir([image_path '\*.jpg']);
file_num=length(files)
names = cell(file_num,1);
F = zeros(file_num,17000, 'single');

for w =1:file_num
    tic
    jpegfilename = [image_path '\' files(w).name];   %% jpegfilename = 'C:\Users\VAZ\Desktop\3.jpg';
    beta = cal_cover_difference_from_QF_3(jpegfilename,QF,num_size);
    f = QF_GFR(jpegfilename, QF, beta);
    F(w,:) = f(:);   %将f以行向量的形式存储， w代表样本的个数。即每个样本的特征通过F对应一行的元素表示 %写成：的形式，就可以用parfor
    names{w} = files(w).name;
    toc
end

save(feature_path,'F','names','-v7.3');
disp('end')