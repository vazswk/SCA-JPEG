% mian_JUNIWARD.m
clear;clc;
%matlabpool 6;

QF = 75;
% QF = 95;
% CAPA = [0.05 0.1 0.2 0.3 0.4];
CAPA = 0.4;
clen = length(CAPA);
Error_JUNIWARD = zeros(1,clen);
Std_JUNIWARD = zeros(1,clen);

indir = 'D:\panyuanfeng\BossBase-10000-jpgQ75';
% feature_cover = 'D:\panyuanfeng\feature\JSRM\Q75\10000_cover_BossBase.mat';
% my_ccJRM_SRM_mex(indir,feature_cover,QF);
feature_cover = 'D:\panyuanfeng\feature\SCA_DCTR\Q75\10000_cover_JUNIWARD_40.mat';
% my_ccJRM(indir,feature_cover,QF);

outdir = 'D:\panyuanfeng\MG_test_stego\Q75\JUNIWARD';
MAP_path = 'D:\panyuanfeng\feature\SCA_DCTR\Q75\Pro_JUNIWARD';

for k = 1:clen
flist = dir([indir '\*.jpg']);
flen = length(flist);
fprintf('%s%d\n', 'the num of the files: ',flen);
rate = single(CAPA(k)); % «∂»Î¬  
out_path = [outdir '\' num2str(rate*100)];
map_path = [MAP_path '\' num2str(rate*100)];
% distortion = zeros(flen,1);
% name = cell(flen,1);

config.STC_h = uint32(0); % 0 for optimal coding simulator otherwise sets STC submatrix height (try 7-12)
config.seed  = int32(0); % random seed

tic
%     parfor i = 1:flen

    for pic_num = 1:10
        beta = cell(flen,1);
        name =  cell(flen,1); 
     for i = (pic_num-1)*1000+1 : pic_num*1000
         
       tic
       fprintf('%d%s\n',i, ['   processing image: ' flist(i).name]);
       in_file  = [indir '\' flist(i).name];
       out_file = [out_path '\' flist(i).name]; 
       map_file = [map_path '\' flist(i).name(1:end-4) '.mat'];
%        dist = J_UNIWARD(in_file, out_file, rate, config);
% %        name{i} = flist(i).name;
% %        distortion(i) = dist;
      [S_STRUCT, CP1, CM1] = J_UNIWARD_matlab(in_file, rate);
      jpeg_write(S_STRUCT, out_file);
      
       beta{i} = single(CP1);  
       name{i} = map_file;      
      toc;
     end
     for j = (pic_num-1)*1000+1 : pic_num*1000                
        Ori_beta = single(beta{j});
        map_file = name{j};
        save(map_file,'Ori_beta');                    
     end
     clear beta;
     clear name;
    end
toc

matlabpool close;
pause(10);
matlabpool 6;
my_SCA_DCTR(indir,feature_cover,map_path,QF);

image_stego = out_path;
feature_path = 'D:\panyuanfeng\feature\SCA_DCTR\Q75\10000_JUNIWARD_sim.mat';
feature_stego = [feature_path(1:end-4) '_' num2str(rate*100) '.mat'];

matlabpool close;
pause(10);
matlabpool 6;
% my_ccJRM(image_stego,feature_stego,QF);
my_SCA_DCTR(image_stego,feature_stego,map_path,QF);
matlabpool close;
pause(10);
[Error_JUNIWARD(k),Std_JUNIWARD(k)] = my_ensemble(feature_cover,feature_stego); 
% delete('F:\panyuanfeng\MG_test_results\*.jpg');
matlabpool 6;
end
% save('Error_JUNIWARD_Sim','Error_JUNIWARD');
% save('Std_JUNIWARD_Sim','Std_JUNIWARD');
%matlabpool close;
Error_JUNIWARD
 