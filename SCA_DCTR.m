function F = SCA_DCTR(IMAGE, QF, Beta)
% -------------------------------------------------------------------------
% Copyright (c) 2016 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Contact:     mboroum1@binghamton.edu   | May 2016
%              tdenema1@binghamton.edu    
%              fridrich@binghamton.edu 
%          
% http://dde.binghamton.edu/download/feature_extractors
% -------------------------------------------------------------------------
% This function extracts selection channel aware DCTR features for 
% steganalysis of JPEG images proposed in [1]. 
% Note: Phil Sallee's Matlab JPEG toolbox(function jpeg_read)is needed for 
% running this function.
%  -------------------------------------------------------------------------
% Input:  IMAGE......path to JPEG image 
%         QF ........JPEG quality factor  
%         Beta ......embedding change probabilities 
% Output: F .........extracted features
% -------------------------------------------------------------------------
% [1] Tomáš Denemark, Mehdi Boroumand and Jessica Fridrich, "Steganalysis
% Features for Content-Adaptive JPEG Steganography", IEEE TRANSACTIONS ON
% INFORMATION FORENSICS AND SECURITY, VOL. 11, NO. 8, AUGUST 2016 
% -------------------------------------------------------------------------

I_STRUCT=jpeg_read(IMAGE);

% number of histogram bins
T = 4;

% compute quantization step based on quality factor
if QF<50,
    q = min(8 * (50 / QF), 100);
else
    q = max(8 * (2 - (QF/50)), 0.2);
end

%computing DCT bases
k=0:7;
l=0:7;
[k,l]=meshgrid(k,l);

A=0.5*cos(((2.*k+1).*l*pi)/16);
A(1,:)=A(1,:)./sqrt(2);
A=A';

DCTbase = cell(8,8);
for mode_r = 1:8
    for mode_c = 1:8
        DCTbase{mode_r,mode_c}  = A(:,mode_r)*A(:,mode_c)';
    end
end

% Compute DCTR locations to be merged
mergedCoordinates = cell(25, 1);
for i=1:5
    for j=1:5
        coordinates = [i,j; i,10-j; 10-i,j; 10-i,10-j];
        coordinates = coordinates(all(coordinates<9, 2), :);
        mergedCoordinates{(i-1)*5 + j} = unique(coordinates, 'rows');
    end
end

% Decompress to spatial domain
QT=I_STRUCT.quant_tables{1};
fun = @(x) idct2(x.data.*QT);
I_spatial = blockproc(I_STRUCT.coef_arrays{1},[8 8],fun);

% Compute features
spatial=zeros(size(Beta));
for k = 1:8
    for l = 1:8
        spatial = spatial + kron(Beta(k:8:end, l:8:end),abs(DCTbase{k,l})*QT(k,l));
    end
end

modeFeaDim = numel(mergedCoordinates)*(T+1);
F = zeros(1, 64*modeFeaDim, 'single');

for mode_r = 1:8
    for mode_c = 1:8
        modeIndex = (mode_r-1)*8 + mode_c;
                  
        % Obtain DCT residual R by convolution between image in spatial domain and the current DCT base
        R = conv2(I_spatial, DCTbase{mode_r,mode_c}, 'valid');
                
        % Quantization, rounding, absolute value, thresholding
        R = abs(round(R / q));      
        R(R > T) = T;
        
        % Selection channel quantity
        quantity = sqrt(conv2(spatial, abs(DCTbase{mode_r, mode_c}), 'valid')); 
        
        % Core of the feature extraction
        for merged_index=1:numel(mergedCoordinates)
            f_merged = zeros(1, T+1, 'single');
            for coord_index = 1:size(mergedCoordinates{merged_index}, 1);
                r_shift = mergedCoordinates{merged_index}(coord_index, 1);
                c_shift = mergedCoordinates{merged_index}(coord_index, 2);
                R_sub = R(r_shift:8:end, c_shift:8:end);
                q_sub = quantity(r_shift:8:end, c_shift:8:end);
                f_merged = f_merged + accumarray(R_sub(:)+1,q_sub(:),[T+1 1])';
            end
            F_index_from = (modeIndex-1)*modeFeaDim+(merged_index-1)*(T+1)+1;
            F_index_to = (modeIndex-1)*modeFeaDim+(merged_index-1)*(T+1)+T+1;
            F(F_index_from:F_index_to) = f_merged / sum(f_merged);
        end
    end
end