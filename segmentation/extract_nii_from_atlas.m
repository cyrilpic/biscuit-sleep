% Extract individual structure layer from AAL3

uiimport('AAL3_1mm.nii.txt')
%%
vol_h = spm_vol('AAL3_1mm.nii');
vol = spm_read_vols(vol_h);
%%
mkdir AAL3_1mm_LR

for i=1:size(AAL31mm, 1)
    num = AAL31mm{i, 1};
    name = char(AAL31mm{i, 2});
    s_vol = vol == num;
    s_vol_h = spm_create_vol(vol_h);
    s_vol_h.fname = ['AAL3_1mm_LR/' name '.nii'];
    s_vol_h.descrip = [s_vol_h.descrip ' ' name];
    spm_write_vol(s_vol_h, uint8(s_vol));
end
%%
mkdir AAL3_1mm

for i=1:size(AAL31mm, 1)
    num = AAL31mm{i, 1};
    name = char(AAL31mm{i, 2});
    if strcmp(name(end-1:end), '_L')
        name = name(1:end-2);
        s_vol = (vol == num) | (vol == num+1);
    elseif strcmp(name(end-1:end), '_R')
        continue
    else
        s_vol = (vol == num);
    end
    s_vol_h = spm_create_vol(vol_h);
    s_vol_h.fname = ['AAL3_1mm/' name '.nii'];
    s_vol_h.descrip = [s_vol_h.descrip ' ' name];
    spm_write_vol(s_vol_h, s_vol);
end