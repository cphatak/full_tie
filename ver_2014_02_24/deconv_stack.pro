@read_dm3
;+
;NAME: DECONV_STACK.pro
;
;USAGE: deconv_stack,fs
;
;PURPOSE: This routine read the stack of through focus series images, and deconvoles each
;         image with the Modulation Transfer Function of the CCD, and saves the stack as 
;         *.decon. The program either rebins the data to 1K by 1K or lets user select a 
;         1K by 1K area from the image.
;
;ARGUMENTS: fs - pointer for file structure generated using rdstack.pro
;
;NEEDS: read_dm3.pro
;
;AUTHOR: CD Phatak, ANL, 09/14/2012.
;
;Modified: Added option of saving the data as raw file for imageJ alignment.
;-
pro deconv_stack, fs, create_raw=create_raw
compile_opt idl2

;get dimensions
dim = fs.dim
;This program works with 1K
;dim = 1024

;close all open files
close,/all

;read iMTF file
print,format="('loading inverse MTF...',$)"
openr,1,'iMTF-ANL-Orius.data'
q=assoc(1,intarr(1))
ndim = (q[0])
ndim = ndim[0]
q=assoc(1,fltarr(ndim,ndim),2)
iMTF = (q[0])
close,1
print,'done'


;read the initial dm3 stack
read_dm3,fs.path+fs.prefix+'.dm3',iminfo
imstack = float(iminfo.image)

;start deconvolution
;the original stack consists of images from -defocus to +defocus
;however for the new deconvolved stack due to legacy coding, the
;stack will be saved as infocus,-defocus, defocus.
deconstack = fltarr(dim,dim,fs.num)

window,0,xsi=dim,ysi=dim,retain=2
for inum = 0,fs.num-1 do begin
  image = imstack[*,*,inum]
  ;check if user wants to rebin the image or select a 1K area.
  if (inum eq 0) then begin
		if (fs.dim ge dim) then begin
		resz = 0
		read,prompt='Do you want to resize the image (1) or select a 1k by 1k area (0)?',resz
		endif
		if (resz eq 0) then begin
		print,'Click on a point around which a 1k by 1k area will be selected'
		sat=0
		while (sat eq 0) do begin
		erase
		tvscl,rebin(imstack[*,*,(fs.num-1)/2],dim,dim)
		cursor,x,y,4,/dev
		plots,x,y,/dev,color=200,psym=7
		scl=fs.dim/dim
			x*=scl & y*=scl
		if (x lt fs.dim/2) then begin
			xs = 0 > (x-dim/2)
			xe = (x+dim/2-1) > (dim-1)
		endif else begin
			xs = (fs.dim-1-dim+1) < (x-dim/2)
			xe = (x+dim/2-1) < (fs.dim-1)
		endelse
		if (y lt fs.dim/2) then begin
			ys = 0 > (y-dim/2)
			ye = (y+dim/2-1) > (dim-1)
		endif else begin
			ys = (fs.dim-1-dim+1) < (y-dim/2)
			ye = (y+dim/2-1) < (fs.dim-1)
		endelse
	    plots,[xs/scl,xe/scl,xe/scl,xs/scl,xs/scl],$
	  		  [ys/scl,ys/scl,ye/scl,ye/scl,ys/scl],/dev,color=200
		read,prompt='Are you satisfied with the selection (0/1):',sat
		cc = tvrd(/order,true=1)
		endwhile
		;save the selected area image
		write_tiff,fs.path+fs.prefix+'_crparea.tiff',cc
		endif
	endif
  
  ;processing the image
  if (inum eq 0) then print,format="('Deconvolving and saving all images ',$)"
  tvscl,rebin(image[*,*],dim,dim)
	xyouts,101,101,'Raw',charsize=2.0,/dev
	wait,0.1
	;removing bad pixels using median filter
	image=median(image,3)
	tvscl,rebin(image[*,*],dim,dim)
	xyouts,101,101,'Bad pixels removed',charsize=2.0,/dev
	empty
	fimage = float(fft(fft(image,-1)*iMTF,1))
	tvscl,rebin(fimage,dim,dim)
	xyouts,101,101,'Deconvolved',charsize=2.0,/dev
	empty
  if (resz eq 1) then begin
		fimage = congrid(fimage,dim,dim,cubic=-0.5)
	endif else begin
		fimage = fimage[xs:xe,ys:ye]
	endelse

  ;store the image
  deconstack[*,*,inum] = fimage
  print,format="('.',$)"
endfor
print,'done.'

;make sure the integrated intensity in each image is the same::
t = total(total(deconstack,1),1)
;print,'image totals : ',t
tm = max(t)
t = tm/t
for kk=0,fs.num-1 do deconstack[0:*,0:*,kk] *= t[kk]


;open the file to store the deconvolved images
openw,10,fs.path+fs.prefix+'.decon'
q=assoc(10,intarr(3))
q[0] = [dim,dim,fs.num]
q=assoc(10,fltarr(dim,dim),12)
;store in-focus image
q[0] = deconstack[*,*,(fs.num-1)/2]
;store the under focus images next
for iunder = 1,(fs.num-1)/2 do q[iunder] = deconstack[*,*,(fs.num-1)/2-iunder]
;store the over focus images
for iover = (fs.num-1)/2+1,fs.num-1 do q[iover] = deconstack[*,*,iover]
close,10

if (keyword_set(create_raw)) then begin
	openw,11,fs.path+fs.prefix+'.raw'
	swap_endian_inplace,deconstack
	writeu,11,deconstack
	close,11
endif

wdelete,0
;update info in the file structure.
scl_fac = dim/fs.dim
fs.dim=dim
if (resz eq 1) then fs.delta *= scl_fac

end
