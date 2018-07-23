files = spm_select([1 Inf],'image','Select con or T images:');
numfiles = size(files,1);

for ii=1:numfiles
    timg = spm_read_vols(spm_vol(files(ii,:))); 
    img(:,:,:,ii)=timg; 
    clear timg;
end

Cimg=zeros(size(img));
Nimg=zeros(size(img));
Pimg=zeros(size(img));
Cimg(~isnan(img))=1;
Nimg(img<0)=1;
Pimg(img>0)=1;
Cimg=sum(Cimg,4);
Pimg=nansum(Pimg,4);
Bimg=zeros(size(Pimg));
hdr=spm_vol('mask.nii');
for ii=1:max(Cimg(:))
    for jj=1:size(img,4)
        ind1=find(Pimg==ii);
        ind2=find(Cimg==jj);
        ind=intersect(ind1,ind2);
        Bimg(ind)=binocdf(ii,jj,.5);
        clear ind1 ind2 ind
    end
end
Bimg=Bimg;
hdr.fname='test3.nii';
spm_write_vol(hdr,Bimg);
clear Bimg
Nimg=nansum(Nimg,4);
Bimg=zeros(size(Nimg));
for ii=1:max(Cimg(:))
    for jj=1:size(img,4)
        ind1=find(Nimg==ii);
        ind2=find(Cimg==jj);
        ind=intersect(ind1,ind2);
        Bimg(ind)=binocdf(ii,jj,.5);
        clear ind1 ind2 ind
    end
end
Bimg=Bimg;
hdr.fname='test3_neg.nii';
spm_write_vol(hdr,Bimg)
