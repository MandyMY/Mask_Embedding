function refocused_Im = refocus(LF, slope,eta)
[U, V, row, col, ~] = size(LF);

M = U * V;
LF_Image = zeros(M, row, col, 3);
k = 0;
for u = 1 : U
    for v = 1 : V
        k = k + 1;
        LF_Image(k, :, :, :) = squeeze(LF(u, v, :, :, :));
 %      imwrite(uint8(squeeze(LF_Image(k,:,:,:))), ['LF_Image', num2str(k), '.png'])
    end
end

shifted_Im = zeros(M, row, col, 3);

UVec =  linspace(-0.5,0.5, U) *(U-1) * slope * eta;
VVec = linspace(-0.5,0.5, V) * (V-1) * slope;
UMat = repmat(UVec', 1, V);
VMat = repmat(VVec, U, 1);
UMat = UMat';
VMat = VMat';
D = [VMat(:), UMat(:)];

for k = 1 : M
    shifted_Im(k,:,:,:) = ImWarp(squeeze(LF_Image(k,:,:,:)), D(k,1), D(k,2));
end

refocused_Im = squeeze(sum(shifted_Im, 1)/M);
