function I = convertSliceToImage(slice,mask_slice)
    I = nan([size(slice),3]);
    hot_idx = slice>0; cold_idx = slice<0;
    hot_rgb = [1,0,0]; cold_rgb = [0,0,1];
    
    for c = 1:3
        tmp = 1-mask_slice;
        tmp(mask_slice) = 0.5;
        tmp(hot_idx) = hot_rgb(c);
        tmp(cold_idx) = cold_rgb(c);
        I(:,:,c) = fliplr(tmp);
    end
    
    I = permute(I,[2,1,3]);
end