;+
;NAME: raw2aligned.pro
;
;CALLING_SEQUENCE: raw2aligned, fs
;
;PURPOSE: This routine converts the raw image stack aligned using ImageJ into
;         a *.aligned stack, to be used by further routines in IDL. The input
;         is the pointer generated using the rdstack.pro routine. The default
;         filename assumed for the raw image stack is tfs.prefix+'_alimJ.raw'
;         and the size is assumed to be 2K by 2K images.
;
;ARGUMENTS: fs - pointer for file structure generated using rdstack.pro
;
;AUTHOR: CD Phatak, ANL, 01/21/2013.
;-
pro raw2aligned, fs
compile_opt idl2

;set the input dimensions to be 2048
idim = fs.dim
id2=idim/2
;initialize the stack variable
imst = fltarr(idim,idim,fs.num)
;open the file for reading
fnraw = fs.prefix+'_alimj.raw'
print,'Reading image stack from raw file '+fnraw
fi = file_info(fs.path+fnraw)
if (fi.exists ne 1) then begin
  print,'The raw file does not exist.'
  goto, cancel
endif
openr,1,fs.path+fnraw
readu,1,imst
close,1
;swap endian since ImageJ uses big endian.
swap_endian_inplace, imst

;generate the mask from the stack
mask = replicate(1.0,idim,idim)
for ii=0,fs.num-1 do begin
  temp = reform(imst[*,*,ii])
  qq = where(temp eq 0,count)
  if (qq[0] ne -1) then begin
    mask[qq]=0.0
    temp[qq]=-1.0
  endif
  imst[*,*,ii]=median(temp,3)
endfor

;resize all the data to 1K by 1K.
dim=2048;fs.dim
d2=dim/2
;imst=imst[id2-d2:id2+d2-1,id2-d2:id2+d2-1,*]
;mask=mask[id2-d2:id2+d2-1,id2-d2:id2+d2-1]
fs.dim=dim
;if (fs.dim ne 1024) then begin
;  fs.dim=1024
;  fs.delta *= 2.0
;endif
;imst = rebin(imst,dim,dim,fs.num)
;mask = rebin(mask,dim,dim)

;show the data
wsi=dim/2
window,0,xsi=wsi,ysi=wsi,title='Aligned data set'
for rpt=0,5 do begin
for ii=0,fs.num-1 do begin 
  tvscl, congrid(imst[*,*,ii]*mask,wsi,wsi)
  wait,0.05
endfor
for ii=fs.num-1,0,-1 do begin
  tvscl,congrid(imst[*,*,ii]*mask,wsi,wsi)
  wait,0.05
endfor
endfor

;store the data into *.aligned file
savfil = fs.path+fs.prefix+'.aligned'
fi = file_info(savfil)
if (fi.exists eq 1) then begin
  ans='y'
  read,prompt='Aligned file already exists. Overwrite? (y/n): ',ans
  if (ans eq 'n') then goto, cancel
endif
print,'Saving the aligned file.'
openw,1,savfil
q=assoc(1,intarr(3))
q[0]=[dim,dim,fs.num]
q=assoc(1,fltarr(dim,dim),12)
;store the in-focus image first
q[0] = imst[*,*,(fs.num-1)/2]
;next store the underfocus images
for ii=1,(fs.num-1)/2 do q[ii] = imst[*,*,(fs.num-1)/2-ii]
;next store the overfocus images
for ii=(fs.num-1)/2+1,fs.num-1 do q[ii]=imst[*,*,ii]
;lastly store the mask
q[fs.num] = float(mask)
close,1

print,'Aligned image stack stored in '+fs.prefix+'.aligned'
wdelete,0

;remove the *.decon files
ans='n'
;check if decon file exists. 
;fi_decon = file_info(fd.path+fd.prefix+'.decon')
;if (fi_decon.exists eq 1) then begin
;read,prompt='Do you want to delete the *.decon file?',ans
;if (ans eq 'y') then spawn,'rm '+fd.path+fd.prefix+'.decon'
;endif

cancel:
end

