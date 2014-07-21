pro MTF_decon,fd,create_raw=create_raw
;
; this routine takes a bunch of 1Kx1K images and deconvolves the 
; Modulation Transfer Function stored in ../MTF/iMTF-ANL-GIF.data
;
; the input to this routine is the output of the rd_fls routine
;
; the output is the file prefix.decon, where prefix is one of the fields
; in fd.
;
;modified: to remove bad pixels due to xray hits, CD, 05/29/07
;modified: to select a 1k by 1k area or resize the entire image to 1k by 1k, CD, 04/01/11

;fd.dim=1024
;dim=1024
dim=fd.dim
close,/all

;resizing
device,cursor_standard=43

; read iMTF file
print,format="('loading inverse MTF ',$)"
;openr,1,'../../mtffiles/iMTF-ANL-Orius.data'
openr,1,'iMTF-ANL-Orius.data'
q=assoc(1,intarr(1))
ndim = (q[0])
ndim = ndim[0]
q=assoc(1,fltarr(ndim,ndim),2)
iMTF = (q[0])
iMTF = rebin(iMTF,dim,dim)
close,1
print,'done'

; deconvolve all images
openw,10,fd.path+fd.prefix+'.decon'
q=assoc(10,intarr(3))
q[0] = [dim,dim]
q=assoc(10,fltarr(dim,dim),12)

;iMTF=rebin(iMTF,dim,dim)
print,format="('Deconvolving and saving all images ',$)"
window,0,xsi=dim,ysi=dim,retain=2

;also correcting the axis angle
axis_angle = 0.0;-31.5
deconstack = fltarr(dim,dim,fd.num)
for i=0L,long(fd.num-1) do begin
	read_dm3,fd.path+fd.names[i],iminfo
	image=float(iminfo.image)
	if (i eq 0) then begin
		if (fd.dim ge dim) then begin
		resz = 0
		read,prompt='Do you want to resize the image (1) or select a 1k by 1k area (0)?',resz
		endif
		if (resz eq 0) then begin
		print,'Click on a point around which a 1k by 1k area will be selected'
		sat=0
		while (sat eq 0) do begin
		erase
		tvscl,rebin(image,dim,dim)
		cursor,x,y,4,/dev
		plots,x,y,/dev,color=200,psym=7
		scl=fd.dim/dim
			x*=scl & y*=scl
		if (x lt fd.dim/2) then begin
			xs = 0 > (x-dim/2)
			xe = (x+dim/2-1) > (dim-1)
		endif else begin
			xs = (fd.dim-1-dim+1) < (x-dim/2)
			xe = (x+dim/2-1) < (fd.dim-1)
		endelse
		if (y lt fd.dim/2) then begin
			ys = 0 > (y-dim/2)
			ye = (y+dim/2-1) > (dim-1)
		endif else begin
			ys = (fd.dim-1-dim+1) < (y-dim/2)
			ye = (y+dim/2-1) < (fd.dim-1)
		endelse
	    plots,[xs/scl,xe/scl,xe/scl,xs/scl,xs/scl],$
	  		  [ys/scl,ys/scl,ye/scl,ye/scl,ys/scl],/dev,color=200
		read,prompt='Are you satisfied with the selection (0/1):',sat
		cc = tvrd(/order,true=1)
		endwhile
		;saave the selected area image
		write_tiff,fd.path+fd.prefix+'_crparea.tiff',cc
		endif
	endif
	;image=congrid(image,dim,dim,cubic=-0.5)
	tvscl,rebin(image[*,*],dim,dim)
	xyouts,101,101,'Raw',charsize=2.0,/dev
	wait,0.1
	;removing bad pixels using median filter
	image=median(image,3)
	tvscl,rebin(image[*,*],dim,dim)
	xyouts,101,101,'Raw: hot spots removed',charsize=2.0,/dev
	empty
	fimage = float(fft(fft(image,-1)*iMTF,1))
	tvscl,rebin(fimage,dim,dim)
	xyouts,101,101,'Deconvolved',charsize=2.0,/dev
	empty
	;shift image intensity such that min. is 0.
	min_img = min(fimage)
	fimage += abs(min_img)
	;correct axis angle
	fimage=rot(fimage,axis_angle,missing=-1.0,cubic=-0.5)
	tvscl,rebin(fimage,dim,dim)
	xyouts,101,101,'Axis corrected.',charsize=2.0,/dev
	if (resz eq 1) then begin
		fimage = congrid(fimage,dim,dim,cubic=-0.5)
	endif else begin
		fimage = fimage[xs:xe,ys:ye]
	endelse
	q[i] = fimage
	deconstack[*,*,i] = fimage
	print,format="('.',$)"
endfor
print,'done'

close,10
print,'deconvolved data stored in '+fd.path+fd.prefix+'.decon'
fd.dim = dim
wdelete,0

if (keyword_set(create_raw)) then begin
	;rearrange the files to be in order for raw file; with zero in the center;
	deconraw_stack = fltarr(dim,dim,fd.num)
	;store the infocus image
	deconraw_stack[*,*,(fd.num-1)/2] = deconstack[*,*,0]
	;next the underfocus images
	for iunder = 1,(fd.num-1)/2 do deconraw_stack[*,*,(fd.num-1)/2-iunder] = deconstack[*,*,iunder]
	;next the overfocus images
	for iover = (fd.num-1)/2+1,fd.num-1 do deconraw_stack[*,*,iover] = deconstack[*,*,iover]
	openw,11,fd.path+fd.prefix+'.raw'
	swap_endian_inplace,deconraw_stack
	writeu,11,deconraw_stack
	close,11
endif


end
