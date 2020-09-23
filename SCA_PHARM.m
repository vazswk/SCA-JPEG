function f = SCA_PHARM(IMAGE, QF, Beta)
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
% This function extracts selection channel aware PHARM features for 
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

settings.T = 2;
settings.nu = 900;
settings.s = 8;
settings.q = (65/4)-(3/20)*QF;

settings.seedIndex = 1;

spatial = DCTtoSpatial( Beta, I_STRUCT.quant_tables{1} );

% Decompress to spatial domain
fun = @(x)x.data .*I_STRUCT.quant_tables{1};
I_spatial = blockproc(I_STRUCT.coef_arrays{1},[8 8],fun);
fun=@(x)idct2(x.data);
X = blockproc(I_spatial,[8 8],fun)+128;

f = struct;
f = all1st(X,1,f);  % 1st order
f = all3rd(X,1,f);  % 3rd order
f = all2x2(X,1,f);

function f = all1st(X,q,f)
    ps = 3; % predictor kernel size
    [M, N] = size(X); [I,J] = deal(2:M-1,2:N-1);
    % Variable names are self-explanatory (R = right, U = up, L = left, D = down)
    [R,U]  = deal(X(I,J+1)-X(I,J),X(I-1,J)-X(I,J));
    f.s1_spam14_R = reshape(ProjHistSpam(R,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
    f.s1_spam14_U = reshape(ProjHistSpam(U,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
end
function f = all3rd(X,q,f)
    ps = 5;
    [M, N] = size(X); [I,J,] = deal(3:M-2,3:N-2);
    [R,U] = deal(-X(I,J+2)+3*X(I,J+1)-3*X(I,J)+X(I,J-1),-X(I-2,J)+3*X(I-1,J)-3*X(I,J)+X(I+1,J));
    f.s3_spam14_R = reshape(ProjHistSpam(R,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
    f.s3_spam14_U = reshape(ProjHistSpam(U,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
end
function f = all2x2(X,q,f)
    ps = 3;
    [M, N] = size(X); [I,J] = deal(2:M-1,2:N-1);
    [Dh, Dv, Dd]  = deal(X(I-1,J-1)+X(I-1,J)-X(I,J-1)-X(I,J), X(I-1,J-1)-X(I-1,J)+X(I,J-1)-X(I,J), -X(I-1,J-1)+X(I-1,J)+X(I,J-1)-X(I,J));
    f.s2x2_spam14_H = reshape(ProjHistSpam(Dh,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
    f.s2x2_spam14_V = reshape(ProjHistSpam(Dv,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
    f.s2x2_spam14_DMaj = reshape(ProjHistSpam(Dd,ps,q),[],1);settings.seedIndex=settings.seedIndex+1;
end

function h = ProjHistSpam(D, kernelSize, centerVal)
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',settings.seedIndex));
    binEdges = [0:settings.T inf];
    binEdges = binEdges * settings.q * centerVal;
    h = zeros(settings.T * settings.nu, 1);
    for projIndex = 1:settings.nu
        Psize = randi(settings.s, 2, 1);
        shift = randi(8, 2, 1);
        
        P = randn(Psize(1), Psize(2));
        n = sqrt(sum(P(:).^2));
        P = P ./ n;
        
        h_proj = zeros(settings.T+1, 1);
                
        % Normal P
        proj = conv2(D, P, 'valid');
        Q = sqrt(conv2(spatial, abs(P), 'valid'));
        Q = Q( 1:size(proj,1), 1:size(proj,2) ); 
        subsamplingRows = shift(1):8:size(proj, 1);
        subsamplingCols = shift(2):8:size(proj, 2);
        proj = proj(subsamplingRows, subsamplingCols);
        Q = Q(subsamplingRows, subsamplingCols);
        h_neigh = max_hist( abs(proj), Q, binEdges );
        h_proj = h_proj + h_neigh';
        
        % Vertical kernel P flip
        proj = conv2(D, flipud(P), 'valid');
        Q = sqrt(conv2(spatial, abs(flipud(P)), 'valid'));
        Q = Q( 1:size(proj,1), 1:size(proj,2) ); 
        subsamplingRows = 1-shift(1)-size(P, 1)+kernelSize:8:size(proj, 1); subsamplingRows = subsamplingRows(subsamplingRows>0);
        subsamplingCols = shift(2):8:size(proj, 2);
        proj = proj(subsamplingRows, subsamplingCols);
        Q = Q(subsamplingRows, subsamplingCols);
        h_neigh = max_hist( abs(proj), Q, binEdges );
        h_proj = h_proj + h_neigh';
        
        % Horizontal kernel P flip
        proj = conv2(D, fliplr(P), 'valid');
        Q = sqrt(conv2(spatial, abs(fliplr(P)), 'valid'));
        Q = Q( 1:size(proj,1), 1:size(proj,2) ); 
        subsamplingRows = shift(1):8:size(proj, 1);
        subsamplingCols = 1-shift(2)-size(P, 2)+kernelSize:8:size(proj, 2); subsamplingCols = subsamplingCols(subsamplingCols>0);
        proj = proj(subsamplingRows, subsamplingCols);
        Q = Q(subsamplingRows, subsamplingCols);
        h_neigh = max_hist( abs(proj), Q, binEdges );
        h_proj = h_proj + h_neigh';
        
        % Horizontal and vertical kernel P flip
        proj = conv2(D, rot90(P,2), 'valid');
        Q = sqrt(conv2(spatial, abs(rot90(P,2)), 'valid'));
        Q = Q( 1:size(proj,1), 1:size(proj,2) ); %determined by the residual filter
        subsamplingRows = 1-shift(1)-size(P, 1)+kernelSize:8:size(proj, 1); subsamplingRows = subsamplingRows(subsamplingRows>0);
        subsamplingCols = 1-shift(2)-size(P, 2)+kernelSize:8:size(proj, 2); subsamplingCols = subsamplingCols(subsamplingCols>0);
        proj = proj(subsamplingRows, subsamplingCols);
        Q = Q(subsamplingRows, subsamplingCols);
        h_neigh = max_hist( abs(proj), Q, binEdges );
        h_proj = h_proj + h_neigh';
        
        h((projIndex-1)*(settings.T+1) + 1:projIndex*(settings.T+1), 1) = h_proj;
    end
end

function h = max_hist( proj, Q, binEdges )
    h = zeros(1,3);
    for i = 1:length(binEdges)-1
        ind =(proj >= binEdges(i) & proj < binEdges(i+1));
        h(i) = sum( Q( ind ) );
    end
end

function spatial = DCTtoSpatial( Beta, QT )
    %computing DCT basis
    k=0:7;
    l=0:7;
    [k,l]=meshgrid(k,l);
    A=0.5*cos(((2.*k+1).*l*pi)/16);
    A(1,:)=A(1,:)./sqrt(2);
    A=A';
    for mode_r = 1:8
        for mode_c = 1:8
            DCTbase{mode_r,mode_c}  = A(:,mode_r)*A(:,mode_c)';
        end
    end
    
    spatial=zeros(size(Beta));
    for k = 1:8
        for l = 1:8
            spatial = spatial + kron(Beta(k:8:end, l:8:end),abs(DCTbase{k,l})*QT(k,l));
        end
    end
end
end
