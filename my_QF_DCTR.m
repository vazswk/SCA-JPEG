function my_QF_DCTR(image_path,feature_path, QF, num_size)
files=dir([image_path '\*.jpg']);
file_num=length(files)
names = cell(file_num,1);
F = zeros(file_num,8000, 'single');

parfor w =1:file_num
    tic
    jpegfilename = [image_path '\' files(w).name];   %% jpegfilename = 'C:\Users\VAZ\Desktop\3.jpg';
    beta = cal_cover_difference_from_QF(jpegfilename,QF,num_size);
    f = SCA_DCTR(jpegfilename, QF, beta);
    F(w,:) = f(:);   %��f������������ʽ�洢�� w���������ĸ�������ÿ������������ͨ��F��Ӧһ�е�Ԫ�ر�ʾ %д�ɣ�����ʽ���Ϳ�����parfor
    names{w} = files(w).name;
    toc
end

save(feature_path,'F','names','-v7.3');
disp('end')