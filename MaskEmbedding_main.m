clear all
clc
idxSave = 0;
path = './LF_original_5/'; %% path='./LF_original_15/' to generate 5*15 LFs with masks
recsize = 896;
FILES = dir(fullfile(path));
FILE_NUM = length(FILES) - 2;

for iScene = 1 : FILE_NUM

    
    UVratio = 1;
    U = 5;
    V = 5;  %% V=15 to generate 5*15 LFs with masks
    LF = [];
    LF_Image = [];
    
    %% Read the input LF and form a tensor
    scene_path = [path, FILES(iScene+2).name, '/'];
    files = dir(fullfile(scene_path,'*.png'));
    for u = 1 : U
        for v = 1 : V
            k = (u-1)*V + v;
            I = im2double(imread([scene_path,'/',files(k).name]));
            [H,W,~] = size(I);
            Smin = min([H,W]);
            down_sample = recsize/Smin + 0.01;
            temp = imresize(I, down_sample);
            [H, W, ~] = size(temp);
            Bh = floor(0.5 * (H - recsize));
            Bw = floor(0.5 * (W - recsize));
            LF(u,v,:,:,:) = temp(Bh:Bh+recsize-1,Bw:Bw+recsize-1,:);
            LF_Image(k, :, :, :) = temp(Bh:Bh+recsize-1,Bw:Bw+recsize-1,:);
        end
    end
    
    [slope_begin, slope_end] = textread([scene_path,'/range.txt'],'%f%f');
    
    %%  Perform Refocusing
    
%             refocus_path = ['./refocus/', num2str(iScene,'%03d')];
%             if exist(refocus_path,'dir')==0
%                 mkdir(refocus_path);
%             end
%             SLOPE = linspace(slope_begin, slope_end, 20);
%             for ns = 1 : length(SLOPE)
%                 slope = SLOPE(ns);
%                 refocused_Im = refocus(LF, slope, UVratio);
%                 imwrite(uint8(255*refocused_Im), [refocus_path, '/', num2str(ns,'%03d'), '.png']);
%             end
    
    %% Add single occlusion
    
    for iter_idx = 1 : 20
        
        LFtemp = LF_Image;
        rx = rand();
        if rx >= 0 && rx < 1/3
            LFtemp(:,:,:,1) = LF_Image(:,:,:,2);
            LFtemp(:,:,:,2) = LF_Image(:,:,:,3);
            LFtemp(:,:,:,3) = LF_Image(:,:,:,1);
        end
        if rx >= 1/3 && rx < 2/3
            LFtemp(:,:,:,1) = LF_Image(:,:,:,3);
            LFtemp(:,:,:,2) = LF_Image(:,:,:,1);
            LFtemp(:,:,:,3) = LF_Image(:,:,:,2);
        end
        LF_Image = LFtemp;
        
        %MaskPath = './MaskRatio/';
        MaskPath = './Mask/';
        maskfiles = dir(fullfile(MaskPath,'*.png'));
        MaskNum = length(maskfiles);
        %iMask = iter_idx;
        iMask = ceil(rand() * MaskNum);  %% Mask selection
        if iMask > MaskNum
            iMask = MaskNum;
        end
        Mask_Im = im2double(imread([MaskPath,'/',maskfiles(iMask).name]));
        [MH, MW, MC] = size(Mask_Im);
        if MC == 1
            Mask_Im = cat(3, Mask_Im, Mask_Im, Mask_Im);
        end
        
        %%% random channel shuffling
        Masktemp = Mask_Im;
        rx = rand();
        if rx >= 0 && rx < 1/3
            Masktemp(:,:,1) = Mask_Im(:,:,2);
            Masktemp(:,:,2) = Mask_Im(:,:,3);
            Masktemp(:,:,3) = Mask_Im(:,:,1);
        end
        if rx >= 1/3 && rx < 2/3
            Masktemp(:,:,1) = Mask_Im(:,:,3);
            Masktemp(:,:,2) = Mask_Im(:,:,1);
            Masktemp(:,:,3) = Mask_Im(:,:,2);
        end
        Mask_Im = Masktemp;
        
        Masksize = max(MH, MW);
        MaskDownRatio = recsize/Masksize - 0.02;
        Mask_Im = imresize(Mask_Im, MaskDownRatio);
        [MH, MW, ~] = size(Mask_Im);
        temp = ones(recsize,recsize,3);
        Bh = ceil(0.5 * (recsize - MH + 1));
        Bw = ceil(0.5 * (recsize - MW + 1));
        temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im;
        Mask_Im = temp;
        %figure; imshow(temp,[])
        
        Mask_bin = double(rgb2gray(Mask_Im) >= 0.9);
        %figure; imshow(Mask_bin,[])
        
        idxSave = idxSave + 1;
        fprintf('Current number of generated light fields is %04d\n', idxSave);
        occ_slope = -10 * rand();
        save_path = ['./LF_generated/', num2str(idxSave,'%04d')];
        if exist([save_path,'/occluded'],'dir')==0
            mkdir([save_path,'/occluded']);
        end
        fid = fopen([save_path,'/occluded/range.txt'],'w');
        fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
        fclose(fid);
        %         fid = fopen(['./LF_original/',num2str(iScene,'%03d'),'/range.txt'],'w');
        %         fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
        %         fclose(fid);
        
        
        UVec =  linspace(-0.5,0.5, U) *(U-1) * occ_slope * UVratio;
        VVec = linspace(-0.5,0.5, V) * (V-1) * occ_slope;
        UMat = repmat(UVec', 1, V);
        VMat = repmat(VVec, U, 1);
        UMat = UMat';
        VMat = VMat';
        D = [VMat(:), UMat(:)];
        Dmax = ceil(max(D(:)));
        
        for k = 1 : U * V
            MaskIm_shifted = ImWarp(Mask_Im, -D(k,1), -D(k,2));
            Mask_bin_shifted = ImWarp(Mask_bin, -D(k,1), -D(k,2));
            LF_Masked = squeeze(LF_Image(k,:,:,:)) .* Mask_bin_shifted + MaskIm_shifted .* (1 - Mask_bin_shifted);
            %Dmax = 0.5*(600-512);
            rectifiedImage = LF_Masked(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
            [H, W, ~] = size(rectifiedImage);
            rectifiedImage = imresize(rectifiedImage, 700/H);
            imwrite(uint8(255*rectifiedImage), [save_path,'/occluded/', num2str(k,'%03d'), '.png']);
            
            if k == ceil(0.5 * U * V)
                GT_Image = squeeze(LF_Image(k,:,:,:));
                rectifiedImage = GT_Image(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                [H, W, ~] = size(rectifiedImage);
                rectifiedImage = imresize(rectifiedImage, 700/H);
                imwrite(uint8(255*rectifiedImage), [save_path,'/groundtruth.png']);
                rectifiedMask = Mask_bin(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                rectifiedMask = imresize(rectifiedMask, 700/H);
                imwrite(uint8(255*rectifiedMask), [save_path,'/mask.png']);
            end
        end
    end
    %end
    %
        %% Add double occlusions
        for iter_idx = 1 : 10
            
            LFtemp = LF_Image;
            rx = rand();
            if rx >= 0 && rx < 1/3
                LFtemp(:,:,:,1) = LF_Image(:,:,:,2);
                LFtemp(:,:,:,2) = LF_Image(:,:,:,3);
                LFtemp(:,:,:,3) = LF_Image(:,:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                LFtemp(:,:,:,1) = LF_Image(:,:,:,3);
                LFtemp(:,:,:,2) = LF_Image(:,:,:,1);
                LFtemp(:,:,:,3) = LF_Image(:,:,:,2);
            end
            LF_Image = LFtemp;
            
            MaskPath = './Mask/';
            maskfiles = dir(fullfile(MaskPath,'*.png'));
            MaskNum = length(maskfiles);
    
            iMask1 = ceil(rand() * MaskNum);  %% Mask selection
            if iMask1 > MaskNum
                iMask1 = MaskNum;
            end
            Mask_Im1 = im2double(imread([MaskPath,'/',maskfiles(iMask1).name]));
            [MH, MW, MC] = size(Mask_Im1);
            if MC == 1
                Mask_Im1 = cat(3, Mask_Im1, Mask_Im1, Mask_Im1);
            end
            Masksize = max(MH, MW);
            MaskDownRatio = recsize/Masksize - 0.02;
            Mask_Im1 = imresize(Mask_Im1, MaskDownRatio);
    
            [MH, MW, ~] = size(Mask_Im1);
            temp = ones(recsize,recsize,3);
            Bh = ceil(0.5 * (recsize - MH + 1));
            Bw = ceil(0.5 * (recsize - MW + 1));
            temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im1;
            Mask_Im1 = temp;
            Mask_bin1 = double(rgb2gray(Mask_Im1) >= 0.9);
    
            iMask2 = ceil(rand() * MaskNum);  %% Mask selection
            if iMask2 > MaskNum
                iMask2 = MaskNum;
            end
            Mask_Im2 = im2double(imread([MaskPath,'/',maskfiles(iMask2).name]));
            [MH, MW, MC] = size(Mask_Im2);
            if MC == 1
                Mask_Im2 = cat(3, Mask_Im2, Mask_Im2, Mask_Im2);
            end
            Masksize = max(MH, MW);
            MaskDownRatio = recsize/Masksize - 0.02;
            Mask_Im2 = imresize(Mask_Im2, MaskDownRatio);
    
            [MH, MW, ~] = size(Mask_Im2);
            temp = ones(recsize,recsize,3);
            Bh = ceil(0.5 * (recsize - MH + 1));
            Bw = ceil(0.5 * (recsize - MW + 1));
            temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im2;
            Mask_Im2 = temp;
            Mask_bin2 = double(rgb2gray(Mask_Im2) >= 0.9);
            
            %%% random channel shuffling
            Masktemp = Mask_Im1;
            rx = rand();
            if rx >= 0 && rx < 1/3
                Masktemp(:,:,1) = Mask_Im1(:,:,2);
                Masktemp(:,:,2) = Mask_Im1(:,:,3);
                Masktemp(:,:,3) = Mask_Im1(:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                Masktemp(:,:,1) = Mask_Im1(:,:,3);
                Masktemp(:,:,2) = Mask_Im1(:,:,1);
                Masktemp(:,:,3) = Mask_Im1(:,:,2);
            end
            Mask_Im1 = Masktemp;
            
            
            Masktemp = Mask_Im2;
            rx = rand();
            if rx >= 0 && rx < 1/3
                Masktemp(:,:,1) = Mask_Im2(:,:,2);
                Masktemp(:,:,2) = Mask_Im2(:,:,3);
                Masktemp(:,:,3) = Mask_Im2(:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                Masktemp(:,:,1) = Mask_Im2(:,:,3);
                Masktemp(:,:,2) = Mask_Im2(:,:,1);
                Masktemp(:,:,3) = Mask_Im2(:,:,2);
            end
            Mask_Im2 = Masktemp;
            
    
            idxSave = idxSave + 1;
            fprintf('Current number of generated light fields is %04d\n', idxSave);
            occ_slope1 = -10 * rand();
            occ_slope2 = -10 * rand() - 10;
            save_path = ['./LF_generated/', num2str(idxSave,'%04d')];
            if exist([save_path,'/occluded'],'dir')==0
                mkdir([save_path,'/occluded']);
            end
            fid = fopen([save_path,'/occluded/range.txt'],'w');
            fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
            fclose(fid);
    %         fid = fopen([scene_path,'/range.txt'],'w');
    %         fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
    %         fclose(fid);
    
            UVec1 =  linspace(-0.5,0.5, U) *(U-1) * occ_slope1 * UVratio;
            VVec1 = linspace(-0.5,0.5, V) * (V-1) * occ_slope1;
            UMat1 = repmat(UVec1', 1, V);
            VMat1 = repmat(VVec1, U, 1);
            UMat1 = UMat1';
            VMat1 = VMat1';
            D1 = [VMat1(:), UMat1(:)];
    
            UVec2 =  linspace(-0.5,0.5, U) *(U-1) * occ_slope2 * UVratio;
            VVec2 = linspace(-0.5,0.5, V) * (V-1) * occ_slope2;
            UMat2 = repmat(UVec2', 1, V);
            VMat2 = repmat(VVec2, U, 1);
            UMat2 = UMat2';
            VMat2 = VMat2';
            D2 = [VMat2(:), UMat2(:)];
            Dmax = max(ceil(max(D1(:))), ceil(max(D2(:))));
    
            for k = 1 : U * V
                MaskIm_shifted1 = ImWarp(Mask_Im1, -D1(k,1), -D1(k,2));
                Mask_bin_shifted1 = ImWarp(Mask_bin1, -D1(k,1), -D1(k,2));
                MaskIm_shifted2 = ImWarp(Mask_Im2, -D2(k,1), -D2(k,2));
                Mask_bin_shifted2 = ImWarp(Mask_bin2, -D2(k,1), -D2(k,2));
    
                LF_Masked = squeeze(LF_Image(k,:,:,:)) .* Mask_bin_shifted1 .* Mask_bin_shifted2 ...
                    + MaskIm_shifted1 .* (1 - Mask_bin_shifted1) .* Mask_bin_shifted2...
                    + MaskIm_shifted2 .* (1 - Mask_bin_shifted2);
    
                %Dmax = 0.5*(600-512);
                rectifiedImage = LF_Masked(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                [H, W, ~] = size(rectifiedImage);
                rectifiedImage = imresize(rectifiedImage, 512/H);
                imwrite(uint8(255*rectifiedImage), [save_path,'/occluded/', num2str(k,'%03d'), '.png']);
    
                if k == ceil(0.5 * U * V)
                    GT_Image = squeeze(LF_Image(k,:,:,:));
                    rectifiedImage = GT_Image(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                    rectifiedImage = imresize(rectifiedImage, 512/H);
                    imwrite(uint8(255*rectifiedImage), [save_path,'/groundtruth.png']);
                    rectifiedMask = Mask_bin1(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:) .* ...
                    Mask_bin2(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                    rectifiedMask = imresize(rectifiedMask, 512/H);
                    imwrite(uint8(255*rectifiedMask), [save_path,'/mask.png']);
                end
            end
        end
    
        %% Add triple occlusions
        for iter_idx = 1 : 10
            
            LFtemp = LF_Image;
            rx = rand();
            if rx >= 0 && rx < 1/3
                LFtemp(:,:,:,1) = LF_Image(:,:,:,2);
                LFtemp(:,:,:,2) = LF_Image(:,:,:,3);
                LFtemp(:,:,:,3) = LF_Image(:,:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                LFtemp(:,:,:,1) = LF_Image(:,:,:,3);
                LFtemp(:,:,:,2) = LF_Image(:,:,:,1);
                LFtemp(:,:,:,3) = LF_Image(:,:,:,2);
            end
            LF_Image = LFtemp;
            
            MaskPath = './Mask/';
            maskfiles = dir(fullfile(MaskPath,'*.png'));
            MaskNum = length(maskfiles);
    
            iMask1 = ceil(rand() * MaskNum);  %% Mask selection
            if iMask1 > MaskNum
                iMask1 = MaskNum;
            end
            Mask_Im1 = im2double(imread([MaskPath,'/',maskfiles(iMask1).name]));
            [MH, MW, MC] = size(Mask_Im1);
            if MC == 1
                Mask_Im1 = cat(3, Mask_Im1, Mask_Im1, Mask_Im1);
            end
            Masksize = max(MH, MW);
            MaskDownRatio = recsize/Masksize - 0.02;
            Mask_Im1 = imresize(Mask_Im1, MaskDownRatio);
    
            [MH, MW, ~] = size(Mask_Im1);
            temp = ones(recsize,recsize,3);
            Bh = ceil(0.5 * (recsize - MH + 1));
            Bw = ceil(0.5 * (recsize - MW + 1));
            temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im1;
            Mask_Im1 = temp;
            Mask_bin1 = double(rgb2gray(Mask_Im1) >= 0.9);
    
            iMask2 = ceil(rand() * MaskNum);  %% Mask selection
            if iMask2 > MaskNum
                iMask2 = MaskNum;
            end
            Mask_Im2 = im2double(imread([MaskPath,'/',maskfiles(iMask2).name]));
            [MH, MW, MC] = size(Mask_Im2);
            if MC == 1
                Mask_Im2 = cat(3, Mask_Im2, Mask_Im2, Mask_Im2);
            end
            Masksize = max(MH, MW);
            MaskDownRatio = recsize/Masksize - 0.02;
            Mask_Im2 = imresize(Mask_Im2, MaskDownRatio);
    
            [MH, MW, ~] = size(Mask_Im2);
            temp = ones(recsize,recsize,3);
            Bh = ceil(0.5 * (recsize - MH + 1));
            Bw = ceil(0.5 * (recsize - MW + 1));
            temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im2;
            Mask_Im2 = temp;
            Mask_bin2 = double(rgb2gray(Mask_Im2) >= 0.9);
    
            iMask3 = ceil(rand() * MaskNum);  %% Mask selection
            if iMask3 > MaskNum
                iMask3 = MaskNum;
            end
            Mask_Im3 = im2double(imread([MaskPath,'/',maskfiles(iMask3).name]));
            [MH, MW, MC] = size(Mask_Im3);
            if MC == 1
                Mask_Im3 = cat(3, Mask_Im3, Mask_Im3, Mask_Im3);
            end
            Masksize = max(MH, MW);
            MaskDownRatio = recsize/Masksize - 0.02;
            Mask_Im3 = imresize(Mask_Im3, MaskDownRatio);
    
            [MH, MW, ~] = size(Mask_Im3);
            temp = ones(recsize,recsize,3);
            Bh = ceil(0.5 * (recsize - MH + 1));
            Bw = ceil(0.5 * (recsize - MW + 1));
            temp(Bh:Bh+MH-1, Bw:Bw+MW-1, :) = Mask_Im3;
            Mask_Im3 = temp;
            Mask_bin3 = double(rgb2gray(Mask_Im3) >= 0.9);
            
            
                        %%% random channel shuffling
            Masktemp = Mask_Im1;
            rx = rand();
            if rx >= 0 && rx < 1/3
                Masktemp(:,:,1) = Mask_Im1(:,:,2);
                Masktemp(:,:,2) = Mask_Im1(:,:,3);
                Masktemp(:,:,3) = Mask_Im1(:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                Masktemp(:,:,1) = Mask_Im1(:,:,3);
                Masktemp(:,:,2) = Mask_Im1(:,:,1);
                Masktemp(:,:,3) = Mask_Im1(:,:,2);
            end
            Mask_Im1 = Masktemp;            
            
            Masktemp = Mask_Im2;
            rx = rand();
            if rx >= 0 && rx < 1/3
                Masktemp(:,:,1) = Mask_Im2(:,:,2);
                Masktemp(:,:,2) = Mask_Im2(:,:,3);
                Masktemp(:,:,3) = Mask_Im2(:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                Masktemp(:,:,1) = Mask_Im2(:,:,3);
                Masktemp(:,:,2) = Mask_Im2(:,:,1);
                Masktemp(:,:,3) = Mask_Im2(:,:,2);
            end
            Mask_Im2 = Masktemp;
            
            Masktemp = Mask_Im3;
            rx = rand();
            if rx >= 0 && rx < 1/3
                Masktemp(:,:,1) = Mask_Im3(:,:,2);
                Masktemp(:,:,2) = Mask_Im3(:,:,3);
                Masktemp(:,:,3) = Mask_Im3(:,:,1);
            end
            if rx >= 1/3 && rx < 2/3
                Masktemp(:,:,1) = Mask_Im3(:,:,3);
                Masktemp(:,:,2) = Mask_Im3(:,:,1);
                Masktemp(:,:,3) = Mask_Im3(:,:,2);
            end
            Mask_Im3 = Masktemp;
            
    
            idxSave = idxSave + 1;
            fprintf('Current number of generated light fields is %04d\n', idxSave);
            occ_slope1 = -10 * rand();
            occ_slope2 = -10 * rand() - 10;
            occ_slope3 = -10 * rand() - 20;
            save_path = ['./LF_generated/', num2str(idxSave,'%04d')];
            if exist([save_path,'/occluded'],'dir')==0
                mkdir([save_path,'/occluded']);
            end
            fid = fopen([save_path,'/occluded/range.txt'],'w');
            fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
            fclose(fid);
    %         fid = fopen([scene_path,'/range.txt'],'w');
    %         fprintf(fid, '%1.1f %1.1f',slope_begin, slope_end);
    %         fclose(fid);
    
            UVec1 =  linspace(-0.5,0.5, U) *(U-1) * occ_slope1 * UVratio;
            VVec1 = linspace(-0.5,0.5, V) * (V-1) * occ_slope1;
            UMat1 = repmat(UVec1', 1, V);
            VMat1 = repmat(VVec1, U, 1);
            UMat1 = UMat1';
            VMat1 = VMat1';
            D1 = [VMat1(:), UMat1(:)];
    
            UVec2 =  linspace(-0.5,0.5, U) *(U-1) * occ_slope2 * UVratio;
            VVec2 = linspace(-0.5,0.5, V) * (V-1) * occ_slope2;
            UMat2 = repmat(UVec2', 1, V);
            VMat2 = repmat(VVec2, U, 1);
            UMat2 = UMat2';
            VMat2 = VMat2';
            D2 = [VMat2(:), UMat2(:)];
    
            UVec3 =  linspace(-0.5,0.5, U) *(U-1) * occ_slope3 * UVratio;
            VVec3 = linspace(-0.5,0.5, V) * (V-1) * occ_slope3;
            UMat3 = repmat(UVec3', 1, V);
            VMat3 = repmat(VVec3, U, 1);
            UMat3 = UMat3';
            VMat3 = VMat3';
            D3 = [VMat3(:), UMat3(:)];
    
            Dmax = max(ceil(max(D2(:))), ceil(max(D3(:))));
    
            for k = 1 : U * V
                MaskIm_shifted1 = ImWarp(Mask_Im1, -D1(k,1), -D1(k,2));
                Mask_bin_shifted1 = ImWarp(Mask_bin1, -D1(k,1), -D1(k,2));
                MaskIm_shifted2 = ImWarp(Mask_Im2, -D2(k,1), -D2(k,2));
                Mask_bin_shifted2 = ImWarp(Mask_bin2, -D2(k,1), -D2(k,2));
                MaskIm_shifted3 = ImWarp(Mask_Im3, -D3(k,1), -D3(k,2));
                Mask_bin_shifted3 = ImWarp(Mask_bin3, -D3(k,1), -D3(k,2));
                LF_Masked = squeeze(LF_Image(k,:,:,:)) .* Mask_bin_shifted1 .* Mask_bin_shifted2 .* Mask_bin_shifted3...
                    + MaskIm_shifted1 .* (1 - Mask_bin_shifted1) .* Mask_bin_shifted2 .* Mask_bin_shifted3 ...
                    + MaskIm_shifted2 .* (1 - Mask_bin_shifted2) .* Mask_bin_shifted3...
                    + MaskIm_shifted3 .* (1 - Mask_bin_shifted3);
    
    
                %Dmax = 0.5*(600-512);
                rectifiedImage = LF_Masked(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                [H, W, ~] = size(rectifiedImage);
                rectifiedImage = imresize(rectifiedImage, 512/H);
                imwrite(uint8(255*rectifiedImage), [save_path,'/occluded/', num2str(k,'%03d'), '.png']);
    
                if k == ceil(0.5 * U * V)
                    GT_Image = squeeze(LF_Image(k,:,:,:));
                    rectifiedImage = GT_Image(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                    rectifiedImage = imresize(rectifiedImage, 512/H);
                    imwrite(uint8(255*rectifiedImage), [save_path,'/groundtruth.png']);
                    rectifiedMask = Mask_bin1(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:)...
                        .* Mask_bin2(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:)...
                        .* Mask_bin3(Dmax+1:end-Dmax,Dmax+1:end-Dmax,:);
                    rectifiedMask = imresize(rectifiedMask, 512/H);
                    imwrite(uint8(255*rectifiedMask), [save_path,'/mask.png']);
                end
            end
        end
    %
end
