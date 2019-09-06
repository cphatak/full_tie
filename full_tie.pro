@joint_entropy
@joint_histogram
@mdg_filter_2d
@interp2d
@read_dm3
@tie_reconsep
@colorwheel
;+
;NAME: full_tie.pro
;
;CALLING SEQUENCE: full_tie
;
;PURPOSE: This program implements a GUI version of the full TIE reconstruction routine including flip/unflip align
;
;
;REQUIRES: joint_entropy.pro, joint_histrogram.pro, mdg_filter_2d.pro, interp2d.pro
;
;VERSION: 3.0 - GUI version.
;
;AUTHOR: CD Phatak, ANL, Jul 21st, 2014.
;
;REVISION NOTES: Sep 8 2014 - Complete check of flip/unlfip align params. Corrected the show img function. Everything seems to be working and fixed until flip/unflip align,
;		              Added the settings control for TIE boxes, linked the do_Tierecon button to tie_Reconsep.pro. this was modified to be a function. the tie_reconsep
;			     was modifieid to work with data supplied by full_tie.pro. NOT TESTED YET!!!.
;-
;
;############################################
FUNCTION getmin_fl,param
;############################################
compile_opt idl2

common fd_data, unflip,flip
common data_common, data

line = findgen(unflip.dim)-float(unflip.dim/2)
bintx=replicate(1.0,unflip.dim)##line
binty=replicate(1.0,unflip.dim)#line

;check params
if (data.nparams eq 4) then begin
	c = rot(flip.al_img,param[2],param[3],missing=0.0,cubic=-0.5)
	c = interp2d(c,bintx,binty,bintx+param[0],binty+param[1],missing=0.0,cubic=-0.5)
endif else begin
	c = rot(flip.al_img,param[2],missing=0.0,cubic=-0.5)
	c = interp2d(c,bintx,binty,bintx*param[3]+binty*param[5]+param[0],bintx*param[6]+binty*param[4]+param[1],cubic=-0.5,missing=0.0)
endelse

;cc = smooth(c,3)
;aa = smooth(a,3)
a = unflip.infimg

;check method
if (data.method eq 0) then begin
	c = mdg_Filter_2d(c,'LoG',sigma=flip.sigma)
	if (data.gtlt eq 1) then c=float(c GT flip.thrs)*flip.commask else c=float(c LT flip.thrs)*flip.commask
	a = mdg_Filter_2d(a,'LoG',sigma=unflip.sigma)
	if (data.gtlt eq 1) then a=float(a GT unflip.thrs)*unflip.commask else a=float(a LT unflip.thrs)*unflip.commask
	;SSD Norm to be minimized
	cm = -total((c-a)^2)
	cm = -total(c*a)
endif else begin
	imgs=bytarr(2,unflip.dim,unflip.dim)
	imgs[0,*,*] = bytscl(a*unflip.commask)
	imgs[1,*,*] = bytscl(c*flip.commask)
	res = joint_entropy(imgs,mutual=cm,/nozero)
endelse
tvscl,rebin(0.5*(a+c)*unflip.commask*flip.commask,unflip.dim/4,unflip.dim/4)
return,-cm
end

;############################################
PRO get_mask_corners_event, ev
;############################################

common data_common,data
common cornercount,ccnt,crnum

  WIDGET_CONTROL, ev.TOP, GET_UVALUE=eventval

  if (ev.type eq 0) then begin 
    if (data.curcorner eq 1) then begin
      data.click[ccnt,0,0] = ev.X
      data.click[ccnt,0,1] = ev.Y
      plots,[0,data.dim-1],[data.click[ccnt,0,1],data.click[ccnt,0,1]],color=200,/dev
      plots,[data.click[ccnt,0,0],data.click[ccnt,0,0]],[0,data.dim-1],color=200,/dev
      data.curcorner=2
    end else begin
      data.click[ccnt,1,0] = ev.X
      data.click[ccnt,1,1] = ev.Y
      plots,[0,data.dim-1],[data.click[ccnt,1,1],data.click[ccnt,1,1]],color=200,/dev
      plots,[data.click[ccnt,1,0],data.click[ccnt,1,0]],[0,data.dim-1],color=200,/dev
      data.curcorner=1
      ccnt += 1
      if (ccnt eq crnum) then WIDGET_CONTROL, ev.TOP, /DESTROY
    endelse
  endif

END


;############################################
pro get_mask_corners,dummy
;############################################

common fd_data, unflip,flip
common data_common,data
common cornercount,ccnt,crnum

  IF(XRegistered("get_mask_corners") NE 0) THEN RETURN
;
  base= WIDGET_BASE(TITLE='click on lower left and upper right mask corners for '+string(crnum,format="(I4)")+' masks', /COLUMN)

  draw = WIDGET_DRAW(base, XSIZE=unflip.dim/4, YSIZE=unflip.dim/4, /BUTTON_EVENTS)

  ccnt = 0
  ; Realize the widget hierarchy.
  WIDGET_CONTROL, base, /REALIZE

  ; Retrieve the widget ID of the draw widget. Note that the widget
  ; hierarchy must be realized before you can retrieve this value.
  WIDGET_CONTROL, draw, GET_VALUE=drawID

  ; Make the draw widget the current IDL drawable area.
  WSET, drawID
  ;setup the flipped image as per the settings
  if (data.revset eq 1) then begin
    imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
    mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
  endif else begin
    imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
    mskb = rot(flip.mask,data.rotang,missing=0.0)
  endelse
  imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
  mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
  im = (unflip.infimg+imb)*unflip.mask*mskb
  tvscl,rebin(im,flip.dim/4,flip.dim/4)
  ; Call XMANAGER to manage the widgets.
  XMANAGER, 'get_mask_corners', base
end

;---------------------------------------------------------------------------------------------------------------------------------------------
;Program that handles the main operations....
pro full_tie_ev, event
;
common fd_data, unflip,flip
common widget_common, widget
common data_common, data
common cornercount,ccnt,crn_num

WIDGET_CONTROL, event.id, GET_UVALUE= eventval ; Get the user value that generated the event...
if (N_ELEMENTS(eventval) eq 0) then return
case eventval of

 'LOADSERIES': begin
               filter = ['*,fls','*.dm3']
		path="/Users/cphatak/ANL_work/"
               res3 = dialog_pickfile(title='Select either *.fls file or DM3 stack of images',get_path=pname,filter=filter,path=path)
		if (res3 eq '') then begin
		 print,'Unable to select the file.'
		 goto, cancel
		endif
		dpos = strpos(res3,'.',/reverse_search)
		suffix = strmid(res3,dpos+1)
		plen = strlen(pname)
		prefix = strlowcase(strmid(res3,plen,dpos-plen))
		problem=0
		;file structure with information about the files.
		ds={num:0,prefix:prefix,path:pname,dim:0L,delta:0.0,fexists:bytarr(20)}
		;Now check if we are loading a fls file or a DM3 stack.
		if (suffix eq 'fls') then begin
			;READING FLS file...
			print,'File name prefix: '+prefix
			openr,1,res3
			num=0 & dim=0
			readf,1,num
			ds.num=num
            defvals=fltarr(num/2)
			print,'No. of images in series:',num
			nm=''
			for i=0,num-1 do begin ;Checking individual images.
			 readf,1,nm
			 fi=file_info(pname+nm)
			 if (fi.exists eq 1) then begin
				ds.fexists[i]=1B
			 	print,'Checking for file ',pname+nm,'-> OK'
			 endif else begin
				ds.fexists[i]=0B
				print,'Checking for file',pname+nm,'does not exist.'
			 endelse
			 ;Get the dimensions and delta
			 if (i eq 0) then begin
			 	read_dm3,pname+nm,iminfo
			 	inf_img=iminfo.image
				delta=iminfo.scale
				if (delta[0] < 0.1) then delta *= 1000.0
				sz=iminfo.dims
			 endif
			endfor
			print,'Image Dimensions are ',sz
			print,'Image pixel resolution (nm/px): ',delta
			ds.dim=sz[0]
			ds.delta=delta[0]
			;Now read all the defocus values
			readf,1,defvals
			;ds.defvals[0:num/2-1]=defvals
			close,1
			if (total(ds.fexists) ne byte(num)) then begin
				print,'One or more files are missing. Please Check'
				problem=1
				goto, cancel
			endif
		endif else begin 
			;READING DM3 STACK.......
			;file structure with information about the files.
			;ds={num:0,prefix:prefix,path:pname,dim:0L,delta:0.0,defvals:fltarr(num/2),fexists:bytarr(num)}
			read_dm3,ds.path+ds.prefix+'.dm3',iminfo
			imstack=iminfo.image
			delta=iminfo.scale
			if (delta[0] < 0.1) then delta *= 1000.0
			sz=size(imstack,/dim)
			if (n_elements(sz) ne 3) then begin
				print,'File does not contain 3D image stack..'
				problem=1
				goto, cancel
			endif
			dimx=sz[0] & dimy=sz[1]
			if (dimx ne dimy) then begin
				print,'Stack does not contain square images..'
				problem=1
				goto, cancel
			endif
			ds.dim=dimx
			ds.num=sz[2]
			ds.delta=delta[0]
			defvals=replicate(delta[2],ds.num/2)
			inf_img = imstack[*,*,ds.num/2]
		endelse
		;Now that the file structure/pointer is created, we put that into the big data structure.
		fd = {filedata, dim:fix(ds.dim), delta:ds.delta, defvals:fltarr(ds.num/2), num:ds.num, $
			prefix:prefix, path:pname, infimg:fltarr(ds.dim,ds.dim), sigma:long(1), $
			thrs:long(0), mask:replicate(1.0,ds.dim,ds.dim),disp_img:fltarr(ds.dim,ds.dim),$
			al_img:fltarr(ds.dim,ds.dim),commask:replicate(1.0,ds.dim,ds.dim),tfsfile:'_tfs.align'} 
		;Check which file was being loaded.
		WIDGET_CONTROL,GET_VALUE=val,widget.imnm
		if (val eq 0) then begin; this is unflip image
			if (problem ne 1) then begin
			 unflip={filedata}
			 unflip.infimg=inf_img
			 data.unfl_exists=byte(1)
			 unflip.disp_img=inf_img
			 unflip.num=ds.num & unflip.dim=ds.dim & unflip.delta=ds.delta &unflip.sigma=1.0
			 unflip.prefix=ds.prefix & unflip.path=ds.path & unflip.defvals=defvals
			 unflip.tfsfile=ds.path+ds.prefix+'_tfs.align'
			 data.dim=unflip.dim
			 WIDGET_CONTROL,SET_VALUE=string(unflip.dim,format='(I4)'),widget.dims
			 WIDGET_CONTROL,SET_VALUE=string(unflip.num,format='(I2)'),widget.numtfs
			 WIDGET_CONTROL,SET_VALUE=string(unflip.delta,format='(F6.2)'),widget.delta
			 WIDGET_CONTROL,SET_VALUE='Unflipped series loaded..',widget.status
			endif
		endif else begin ;this is flip image
			if (problem ne 1) then begin
			 flip={filedata}
			 flip.infimg=inf_img
			 flip.disp_img=inf_img
			 data.fl_exists=byte(1)
			 flip.num=ds.num & flip.dim=ds.dim & flip.delta=ds.delta &flip.sigma=1.0
			 flip.prefix=ds.prefix & flip.path=ds.path & flip.defvals=defvals
			 flip.tfsfile=ds.path+ds.prefix+'_tfs.align'
			 data.flipflapfile=flip.path+flip.prefix+'_flipunflip.align'
			 data.dim=flip.dim
			 WIDGET_CONTROL,SET_VALUE=string(flip.dim,format='(I4)'),widget.dims
			 WIDGET_CONTROL,SET_VALUE=string(flip.num,format='(I2)'),widget.numtfs
			 WIDGET_CONTROL,SET_VALUE=string(flip.delta,format='(F6.2)'),widget.delta
			 WIDGET_CONTROL,SET_VALUE='Flipped series loaded..',widget.status
			 endif
		endelse
		;Update the values to display and show the infocus image
		wset,widget.drawid
		tvscl,rebin(inf_img,widget.wsi,widget.wsi)
                range='['+string(min(inf_img))+'-'+string(max(inf_img))+']'
		WIDGET_CONTROL,SET_VALUE=string(0),widget.defstep
		WIDGET_CONTROL,SET_VALUE=string(0),widget.showimgnum
		WIDGET_CONTROL,SET_VALUE=range,widget.range
		WIDGET_CONTROL,SENSITIVE=1,widget.raw2align
		WIDGET_CONTROL,SENSITIVE=1,widget.decon
		if (data.fl_exists+data.unfl_exists eq byte(2)) then begin
			WIDGET_CONTROL,SENSITIVE=1,widget.dotfs
			WIDGET_CONTROL,SENSITIVE=1,widget.doseries
		endif
		cancel:
		if (problem eq 1) then WIDGET_CONTROL,SET_VALUE='Failed to load files..',widget.status
		endcase
 'RAW2ALIGN' : begin
	       ;First check which image series we are dealing with
		WIDGET_CONTROL,GET_VALUE=val,widget.imnm
		if (val eq 0) then fd=unflip else fd=flip
		fnraw = fd.prefix+'_alimj.raw'
		fi=file_info(fd.path+fnraw)
		if (fi.exists ne 1) then begin
			WIDGET_CONTROL,SET_VALUE='Raw image stack not found..',widget.status
			goto, cancel2
		endif
		WIDGET_CONTROL,SET_VALUE='Reading the aligned raw image stack..',widget.status
		imst = fltarr(fd.dim,fd.dim,fd.num)
		;Read the file
		openr,1,fd.path+fnraw
		readu,1,imst
		close,1
		;swap endianness since ImageJ uses big endian
		swap_endian_inplace, imst
		;Generate the mask from the stack.
		mask=replicate(1.0,fd.dim,fd.dim)
		for ii=0,fd.num-1 do begin
			temp=reform(imst[*,*,ii])
			qq=where(temp eq 0,count)
			if (qq[0] ne -1) then begin
				mask[qq]=0.0
				temp[qq]=-1.0
			endif
			imst[*,*,ii] = temp
		endfor
		;show the aligned stack.
		wset,widget.drawID
		for rpt=0,5 do begin
			for ii=0,fd.num-1 do begin
				tvscl,congrid(reform(imst[*,*,ii])*mask,widget.wsi,widget.wsi)
				wait,0.05
			endfor
			for ii=fd.num-1,0,-1 do begin
				tvscl,congrid(reform(imst[*,*,ii])*mask,widget.wsi,widget.wsi)
				wait,0.05
			endfor
		endfor
		;Saving the data...
		savfil = fd.path+fd.prefix+'.aligned'
		fi=file_info(savfil)
		if (fi.exists eq 1) then begin
			ans='y'
			read,prompt='Aligned file already exists. Overwrite (y/n):',ans
			if (ans eq 'n') then goto, cancel2
		endif
		print,'Saving the aligned file...'
		openw,1,savfil
		q=assoc(1,intarr(3))
		q[0]=[fd.dim,fd.dim,fd.num]
		q=assoc(1,fltarr(fd.dim,fd.dim),12)
		;store the in-focus image first
		q[0] = imst[*,*,(fd.num-1)/2]
        ;update the infocus image in the pointer data
        if (val eq 0) then unflip.infimg = imst[*,*,(fd.num-1)/2] else flip.infimg = imst[*,*,(fd.num-1)/2]
		;next store the underfocus images
		for ii=1,(fd.num-1)/2 do q[ii] = imst[*,*,(fd.num-1)/2-ii]
		;next store the overfocus images
		for ii=(fd.num-1)/2+1,fd.num-1 do q[ii]=imst[*,*,ii]
		;lastly store the mask
		q[fd.num] = float(mask)
		close,1
		if (val eq 0) then unflip.mask=mask else flip.mask=mask
		WIDGET_CONTROL,SET_VALUE='Raw image stack converted and stored as'+fd.prefix+'.aligned',widget.status
		cancel2:
		endcase
 'DECONV'    : begin
		;Check which image series
		WIDGET_CONTROL,GET_VALUE=val,widget.imnm
		if (val eq 0) then fs=unflip else fs=flip
		close,/all
		;read the iMTF file
		WIDGET_CONTROL,SET_VALUE='Loading inverse MTF....',widget.status
		openr,1,'iMTF-ANL-Orius.data'
		q=assoc(1,intarr(1))
		ndim=q[0]
		ndim=ndim[0]
		q=assoc(1,fltarr(ndim,ndim),2)
		iMTF=q[0]
		close,1
		;if (fs.dim gt 2048) then iMTF=1.0
		iMTF=1.0
		;open the aligned stack of images
		fi=file_info(fs.path+fs.prefix+'.aligned')
		if (fi.exists ne 1) then begin
			WIDGET_CONTROL,SET_VALUE='Aligned series not found. Run TFS alignment first.',widget.status
			goto, cancel
		endif
		fi=file_info(fs.path+fs.prefix+'.decon')
		if (fi.exists eq 1) then begin
			ans='y'
			read,prompt='Deconvolved Image stack exits. Overwrite (y/n):',ans
			if (ans eq 'n') then goto, cancel3
		endif
		openr,11,fs.path+fs.prefix+'.aligned'
		q11=assoc(11,fltarr(fs.dim,fs.dim),12)
		openw,12,fs.path+fs.prefix+'.decon'
		q12=assoc(12,intarr(3))
		q12=[fs.dim,fs.dim,fs.num]
		q12=assoc(12,fltarr(fs.dim,fs.dim),12)
		WIDGET_CONTROL,SET_VALUE='Deconvolving images and saving...',widget.status
		wset,widget.drawID
		for ii=0,fs.num-1 do begin
			im=q11[ii] ;from the aligned stack.
			tvscl,congrid(im,widget.wsi,widget.wsi)
			xyouts,101,101,'RAW',charsize=2.0,/dev
			wait,0.1
			;Removing bas pixels/zingers using median filter.
			im=median(im,3)
			tvscl,congrid(im,widget.wsi,widget.wsi)
			xyouts,101,101,'RAW: Hot spots removed',charsize=2.0,/dev
			fim=float(fft(fft(im,-1)*iMTF,1))
			tvscl,congrid(fim,widget.wsi,widget.wsi)
			xyouts,101,101,'Deconvolved',charsize=2.0,/dev
			q12[ii]=fim
			print,format="('.',$)"
            ;update the infocus image in the pointer data
            if (ii eq 0) then begin
                if (val eq 0) then unflip.infimg = fim else flip.infimg = fim
            endif
		endfor
		print,'Done'
		q12[fs.num]=q11[fs.num] ;save the mask as well..
		;close files
		close,11
		close,12
		cancel3:
		endcase
 'SHOWIMGNM' : begin
	       ;Check which series
		WIDGET_CONTROL,GET_VALUE=val,widget.imnm
		WIDGET_CONTROL,GET_VALUE=val2,widget.showimgnum
		num_img=fix(val2[0])	
	  	if (val eq 0) then fs=unflip else fs=flip
		fname = fs.path+fs.prefix+'.decon'
		fi=file_info(fname)
		if (fi.exists ne 1) then begin
			fname = fs.path+fs.prefix+'.aligned'
			fi=file_info(fname)
			if (fi.exists ne 1) then begin
				WIDGET_CONTROL,SET_VALUE='Run raw2aligned before proceeding..',widget.status
				goto, cancel4
			endif
		endif
		openr,1,fname
		q=assoc(1,fltarr(fs.dim,fs.dim),12)
		fs.disp_img = q[num_img]
		tempmask=q[fs.num]
		close,1
		;update the image in the original file data pointer
		if (val eq 0) then unflip.disp_img = fs.disp_img else flip.disp_img=fs.disp_img
		;Update the image, sigma, range and defocus value.
		WIDGET_CONTROL,SET_VALUE=string(fs.sigma,format='(F6.2)'),widget.sigma
		imf = mdg_filter_2d(fs.disp_img,'LoG',sigma=float(fs.sigma))
		range='['+string(min(imf))+'-'+string(max(imf))+']'
		WIDGET_CONTROL,SET_VALUE=range,widget.range
		wset,widget.drawID
		tvscl,congrid(fs.disp_img,widget.wsi,widget.wsi)
		if (num_img eq 0) then defval=0.0 else begin
			if (num_img gt fs.num/2) then defval=fs.defvals[num_img-1-fs.num/2] else $
				defval=(-1.0)*fs.defvals[num_img-1]
		endelse
		WIDGET_CONTROL,SET_VALUE=string(defval,format='(F10.2)'),widget.defstep
		cancel4:
		endcase
 'SERIESTYPE' : begin
		;Check which series
		WIDGET_CONTROL,GET_VALUE=val,widget.imnm
		if (val eq 0) then fs=unflip else fs=flip
		;Get the value selected
		WIDGET_CONTROL,GET_VALUE=sval,widget.stype
		case sval of
			0: defvals=fs.defvals
			1: defvals=(findgen(fs.num/2)+1)*fs.defvals[0]
			2: defvals = (findgen(fs.num/2)+1)^2*fs.defvals[0]
		endcase
		WIDGET_CONTROL,SET_VALUE='Defocus values: '+string(defvals),widget.status
		if (val eq 0) then unflip.defvals=defvals else flip.defvals=defvals
		endcase
 'NUMMASK'   : begin 
               WIDGET_CONTROL, GET_VALUE=val, widget.nummask
               data.nummask = fix(val(0))
               WIDGET_CONTROL, SET_VALUE=string(val,format='(I3)'), widget.nummask
               endcase
  'IMNM'     : begin
               WIDGET_CONTROL,GET_VALUE=val, widget.imnm
               if (val eq 0) then begin   ;This is unflip image
                 wset,widget.drawid
		 if (data.unfl_exists eq byte(1)) then begin
                 tvscl,rebin(unflip.infimg,widget.wsi,widget.wsi)
                 WIDGET_CONTROL,SET_VALUE=string(unflip.sigma,format='(F6.2)'),widget.sigma
                 WIDGET_CONTROL,SET_VALUE=string(unflip.thrs,format='(F6.2)'),widget.thrs
                 im=unflip.infimg
                 im=mdg_filter_2d(im,'LoG',sigma=unflip.sigma)
                 range='['+string(min(im))+'-'+string(max(im))+']'
                 WIDGET_CONTROL,SET_VALUE=range,widget.range
		 endif else WIDGET_CONTROL,SET_VALUE='Unflipped series not loaded..',widget.status
               endif else begin     ;Then it is flip image
                 wset,widget.drawid
		 if (data.fl_exists eq byte(1)) then begin
                 tvscl,rebin(flip.infimg,widget.wsi,widget.wsi)
                 WIDGET_CONTROL,SET_VALUE=string(flip.sigma,format='(F6.2)'),widget.sigma
                 WIDGET_CONTROL,SET_VALUE=string(flip.thrs,format='(F6.2)'),widget.thrs
                 im=flip.infimg
                 im=mdg_filter_2d(im,'LoG',sigma=flip.sigma)
                 range='['+string(min(im))+'-'+string(max(im))+']'
                 WIDGET_CONTROL,SET_VALUE=range,widget.range
		 endif else WIDGET_CONTROL,SET_VALUE='Flipped series not loaded..',widget.status
               endelse
               endcase
  'SIGMA'   :  begin
               WIDGET_CONTROL,GET_VALUE=val,widget.sigma
               ;Check which sigma is being updated
               WIDGET_CONTROL,GET_VALUE=val2,widget.imnm
               if (val2 eq 0) then unflip.sigma=float(val[0]) else flip.sigma=float(val[0])
               WIDGET_CONTROL,SET_VALUE=string(val,format='(F6.2)'),widget.sigma
               ;update the intensity range after filtering
               if (val2 eq 0) then fs=unflip else fs=flip
               ;get the image number being displayed
	       im=fs.disp_img
	       imf = mdg_filter_2d(im,'LoG',sigma=float(val[0]))
               range='['+string(min(imf))+'-'+string(max(imf))+']'
               WIDGET_CONTROL,SET_VALUE=range,widget.range
               endcase
  'THRS'    :  begin
               WIDGET_CONTROL,GET_VALUE=val,widget.thrs
               ;Check which sigma is being updated
               WIDGET_CONTROL,GET_VALUE=val2,widget.imnm
               if (val2 eq 0) then unflip.thrs=float(val[0]) else flip.thrs=float(val[0])
               WIDGET_CONTROL,SET_VALUE=string(val,format='(F6.2)'),widget.thrs
               endcase
'DOGDISPLAYGT' : begin
                ;Check which image series is being displayed
                WIDGET_CONTROL,GET_VALUE=val,widget.imnm
                if (val eq 0) then fs=unflip else fs=flip
	        im=fs.disp_img
                  imf = mdg_filter_2d(im,'LoG',sigma=fs.sigma)*fs.mask
		  data.gtlt = 1
                  wset,widget.drawid
                  tvscl,rebin(imf GT fs.thrs,widget.wsi,widget.wsi)
		endcase
'DOGDISPLAYLT' : begin
                ;Check which image is being displayed
                WIDGET_CONTROL,GET_VALUE=val,widget.imnm
                if (val eq 0) then fs=unflip else fs=flip
	        im=fs.disp_img
                  imf = mdg_filter_2d(im,'LoG',sigma=fs.sigma)*fs.mask
                  data.gtlt = -1
                  wset,widget.drawid
                  tvscl,rebin(imf LT fs.thrs,widget.wsi,widget.wsi)
		endcase
 'NPARAMS'    : begin
                WIDGET_CONTROL, GET_VALUE=val,widget.nparams
                if (val eq 0) then data.nparams=4 else data.nparams=7
                endcase
 'METHOD'     : begin
	              WIDGET_CONTROL,GET_VALUE=val,widget.method
	              if (val eq 0) then data.method=0 else data.method=1
	              endcase
 'SHX'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.shx
                data.shx=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                wset,widget.drawid
                tvscl,congrid(unflip.infimg+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.shx
                endcase
 'SHY'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.shy
                data.shy=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                wset,widget.drawid
                tvscl,congrid(unflip.infimg+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.shy
                endcase
 'ROTANG'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.rotang
                data.rotang=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                wset,widget.drawid
                tvscl,rebin(unflip.infimg+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.rotang
                endcase
 'REVSET'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.revset
                if (val eq 0) then data.revset=0 else data.revset=1
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                wset,widget.drawid
                tvscl,congrid(unflip.infimg+imb,widget.wsi,widget.wsi)
                endcase
 'SHOWSUM'        : begin
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                wset,widget.drawid
                tvscl,rebin(unflip.infimg+imb,widget.wsi,widget.wsi)
                endcase
 'MKMASK'     : begin
		 ;Check which image series is being displayed
                WIDGET_CONTROL,GET_VALUE=val,widget.imnm
                if (val eq 0) then fd=unflip else fd=flip

                crn_num = data.nummask
	              get_mask_corners

	              ; create the mask 
	              wexp = 200.D0
	              q = fltarr(fd.dim,fd.dim)
	              for ii=0,data.nummask-1 do q[4*data.click[ii,0,0]:4*data.click[ii,1,0],4*data.click[ii,0,1]:4*data.click[ii,1,1]] = 1.0

	              ; define the coordinate arrays used to compute the Gaussian
	              x = findgen(fd.dim)-float(fd.dim/2)
	              y = findgen(fd.dim)-float(fd.dim/2)
	     
	              ; compute the 2D Gaussian and its integral (to normalize afterwards)
	              !EXCEPT=0
	              s = exp(-x^2/wexp) # exp(-y^2/wexp)
	              norm = total(s)
	     
	              ; FFT both arrays
	              fq = fft(q,-1,/double)
	              fs = fft(s,-1,/double)
	     
	              ; multiply both  (to compute their convolution) and correct amplitude; then shift to center
	              fd.commask = float(shift(float(fft(fq*fs,1,/double)),fd.dim/2,fd.dim/2) * fd.dim * fd.dim / norm)
	              !EXCEPT=1
	              wset,widget.drawID
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(flip.mask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(flip.mask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,/grid)
                fd.commask *= (unflip.mask*mskb)
	              im = (unflip.infimg+imb)*fd.commask
		flip.commask=fd.commask
		unflip.commask=flip.commask
	              tvscl,rebin(im,widget.wsi,widget.wsi)
	              ; once the filter is defined, we can proceed with the alignment
	              WIDGET_CONTROL, SENSITIVE=1, widget.doseries
	              endcase
  'DOALIGN'   : begin
		fi=file_info(data.flipflapfile)
		ans='y'
		if (fi.exists eq 1) then read,ans,prompt='Flip/Unflip align file exists. Overwrite? (y/n):'
		if (ans eq 'y') then begin
              ;Now to begin the actual alignment process
                ;Unflipped image
                aa = smooth(unflip.infimg,3)
                aa = float(mdg_Filter_2d(aa,'LoG',sigma=unflip.sigma))
		if (data.gtlt eq 1) then aa = aa GT unflip.thrs else aa = aa LT unflip.thrs
		tvscl,rebin(aa,widget.wsi,widget.wsi)
		fa = fft(aa*unflip.commask,-1)
		
                ;Flipped image
                ;correct the orientation of the flipped image
                if (data.revset eq 1) then begin
                  imb=reverse(rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0),2)
                endif else begin
                  imb = rot(flip.infimg,data.rotang,cubic=-0.5,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                flip.al_img = imb
                imb = smooth(imb,3)
                imb = float(mdg_filter_2d(imb,'LoG',sigma=flip.sigma))
                if (data.gtlt eq 1) then imb = imb GT flip.thrs else imb = imb LT flip.thrs
                tvscl,rebin(imb,widget.wsi,widget.wsi)
		fb = conj(fft(imb*flip.commask,-1))
		ab = smooth(abs(fft(fa*fb,1)),5)
		tvscl,rebin(shift(alog(ab+0.01),unflip.dim/2,unflip.dim/2),widget.wsi,widget.wsi)
		;determine the max
		qmax=where(ab eq max(ab))
      		roughx = qmax[0] mod unflip.dim
        	roughy = qmax[0]/unflip.dim
        	if (roughx gt unflip.dim/2) then roughx = roughx - unflip.dim
        	if (roughy gt unflip.dim/2) then roughy = roughy - unflip.dim
        	plots,(unflip.dim/2+[roughx-3,roughx+3])/2,(unflip.dim/2+[roughy,roughy])/2,/dev,color=230
        	plots,(unflip.dim/2+[roughx,roughx])/2,(unflip.dim/2+[roughy-3,roughy+3])/2,/dev,color=230
		;start sub-pixel
		ftol=1.0e-4
		stat_updt = 'Rough shifts= '+string(-roughx)+string(-roughy)
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		if (data.nparams eq 4) then begin
			pa = [-roughx,-roughy,0.0,1.0]
			scl = [2.0,2.0,2.0,0.2]
		endif else begin
			pa = [-roughx,-roughy,0.0,1.0,1.0,0.0,0.0]
			scl = [2.0,2.0,2.0,0.1,0.1,0.1,0.1]
		endelse
		iter=0
		itmax=3000
		fmin=amoeba(ftol,function_name='getmin_fl',function_value=fv,ncalls=iter,nmax=itmax,p0=pa,scale=scl)
		print,'Iterations=',iter
		print,'Fits=',fmin
		;resetting the file
		openw,1,data.flipflapfile
		close,1
		openw,2,data.flipflapfile,/append
		printf,2,data.nparams
		printf,2,[data.shx,data.shy,data.rotang,data.revset]
		printf,2,fmin
		close,2
		;update the params in data structure
		if (data.nparams eq 4) then data.params[0:3]=fmin else data.params[0:6]=fmin
		stat_updt = 'Alignment complete. Params saved in '+data.flipflapfile+'. Press Show Align to see the alignment.'
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		endif else begin;end of ans='y'
		;read the alignment data from the previous file....
		openr,1,data.flipflapfile
		tmp=0
		readf,1,tmp
		if (tmp eq 4) then params=fltarr(4) else params=fltarr(7)
		al_set=fltarr(4)
		readf,1,al_set
		readf,1,params
		close,1
		print,al_set
		print,params
		data.shx=al_Set[0]
		data.shy=al_set[1]
		data.rotang=al_set[2]
		data.revset=al_set[3]
		if (tmp eq 4) then data.params[0:3]=params else data.params[0:6]=params
		WIDGET_CONTROL,SET_VALUE=string(al_set[0],format='(I5)'),widget.shx
 		WIDGET_CONTROL,SET_VALUE=string(al_set[1],format='(I5)'),widget.shy
  		WIDGET_CONTROL,SET_VALUE=string(al_set[2],format='(I5)'),widget.rotang
  		WIDGET_CONTROL,SET_VALUE=fix(al_set[3]),widget.revset
		WIDGET_CONTROL,SET_VALUE='Alignment file exists. Pre-loaded the values. Press Show Align to check.',widget.status
		endelse	
		;Apply the alignment and save the image into data structure
		line = findgen(flip.dim)-float(flip.dim/2)
		bintx = replicate(1.0,flip.dim)##line
		binty = replicate(1.0,flip.dim)#line
		if (data.revset eq 1) then begin
			imb = reverse(rot(flip.infimg, data.rotang, cubic=-0.5,missing=0.0),2)
		endif else begin
			imb = rot(flip.infimg, data.rotang, cubic=-0.5, missing=0.0)
		endelse
		imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
		if (data.nparams eq 4) then begin
			imb = rot(imb,data.params[2],data.params[3],cubic=-0.5,missing=0.0)
			imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+params[1],cubic=-0.5,missing=0.0)
		endif else begin
			imb = rot(imb,data.params[2],cubic=-0.5,missing=0.0)
			imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
						bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=0.0,cubic=-0.5)
		endelse
		flip.al_img = imb
        print,data.params
		WIDGET_CONTROL,SENSITIVE=1,widget.showseries
		WIDGET_CONTROL,SENSITIVE=1,widget.goalign
		
		endcase
  'SHOWALIGN' : begin
		wset,widget.drawid
		for cnt = 0,10 do begin
			tvscl,rebin(unflip.infimg*unflip.commask,widget.wsi,widget.wsi)
			wait,0.1
			tvscl,rebin(flip.al_img*flip.commask,widget.wsi,widget.wsi)
			wait,0.1
		endfor
		endcase
  'GOALIGN'   : begin
		stat_updt = 'Applying Alignment to Flipped Images...'
		print,data.params
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		;Input data
		openr,10,flip.path+flip.prefix+'.decon'
		numim = flip.num
		q10 = assoc(10,fltarr(flip.dim,flip.dim),12)
        dec_msk = q10[flip.num]
		;Output data
		openw,20,flip.path+flip.prefix+'_final.aligned'
		q20=assoc(20,intarr(3))
		q20[0]=[flip.dim,flip.dim,numim]
		q20=assoc(20,fltarr(flip.dim,flip.dim),12)
		;Mask from flipped alignment
		mask = replicate(1.0,flip.dim,flip.dim)*flip.mask*dec_msk
		line = findgen(flip.dim)-float(flip.dim/2)
		bintx = replicate(1.0,flip.dim)##line
		binty = replicate(1.0,flip.dim)#line
		for ii=0,numim-1 do begin
			print,'Image #',ii
			imb = q10[ii]
			;Correct the flip...
			 if (data.revset eq 1) then begin
                  		imb=reverse(rot(imb,data.rotang,cubic=-0.5,missing=-1.0),2)
                	 endif else begin
                  		imb = rot(imb,data.rotang,cubic=-0.5,missing=-1.0)
                	 endelse
                	 imb = interpolate(imb,findgen(flip.dim)-data.shx,findgen(flip.dim)-data.shy,missing=-1.0,cubic=-0.5,/grid)
			 ;Apply alignment
			 if (data.nparams eq 4) then begin
				imb = rot(imb,data.params[2],data.params[3],cubic=-0.5,missing=-1.0)
				imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+data.params[1],cubic=-0.5,missing=-1.0)
		         endif else begin
				imb = rot(imb,data.params[2],cubic=-0.5,missing=-1.0)
				imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
						bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=-1.0,cubic=-0.5)
		 	 endelse
			 ;modify the mask
			 qq=where(imb lt 0.0)
			 if (qq[0] ne -1) then mask[qq]=0.0
			 ;store the image
			 q20[ii]=imb
			 ;Show the image
			 wset,widget.drawid
			 tvscl,rebin(imb,widget.wsi,widget.wsi)
		endfor
		print,'----> Done..'
		;save mask
		q20[numim]=mask
        tvscl,rebin(mask,widget.wsi,widget.wsi)
		;close files
		close,10
		close,20
        WIDGET_CONTROL,SENSITIVE=1,widget.dotierecon
		endcase
  'EE'	      : begin
		WIDGET_CONTROL,get_value=val,widget.ee
		data.ee=float(val[0])
		WIDGET_CONTROL,set_Value=string(val,format='(E10.2)'),widget.ee
		endcase
  'CS'	      : begin
		WIDGET_CONTROL,get_Value=val,widget.cs
		data.cs=float(val[0])
		WIDGET_CONTROL,set_Value=string(val,format='(E10.2)'),widget.cs
		endcase
 'TIESETTING' : begin
		WIDGET_CONTROL,get_Value=val,widget.tiesetting
		data.tiesetting=fix(val)
        print,data.tiesetting
		endcase
 'TIKVAL'     : begin
		WIDGET_CONTROL,get_Value=val,widget.tikval
		data.tikval=float(val[0])
		WIDGET_CONTROL,set_Value=string(val,format='(F6.2)'),widget.tikval
		endcase
 'METHODS'    : begin
		WIDGET_CONTROL,get_Value=val,widget.tiemethod
		data.tiemethod=fix(Val)
        print,data.tiemethod
		endcase
 'DERIVMETHOD': begin
		WIDGET_CONTROL,get_value=val,widget.deriv
		data.deriv=fix(Val)
		endcase
 'DOTIERECON' : begin
		;Compile and output all the parameters that will be used for performing TIE_reconsep
		print,'-------Parameters for TIE_RECON-------------'
		print,'Image dimensions: '+string(data.dim,format='(I4)')
		print,'Pixel size: '+string(unflip.delta,format='(F5.2)')
		print,'Acc. Voltage (V): '+string(data.ee,format='(E10.2)')
		print,'Sph. Aberration (nm): '+string(data.cs,format='(E10.2)')
        if data.tiemethod eq 0 then met='Laplacian' else met='Inverse Gradient'
        print,'Method: '+met
        print,'TIE setting: ',data.tiesetting
        print,'Tikhonov value: '+string(data.tikval,format='(I02)')
        print,'------------------------------------------------'
		WIDGET_CONTROL,set_Value='Calling tie_reconsep.pro to perform the phase reconstruction.',widget.status
		result=tie_reconsep(unflip=unflip,flip=flip,data=data)
		if (result eq 1) then stat = 'Phase reconstruction complete without errors.' else stat='Error in Phase reconstruction..' 
		WIDGET_CONTROL,set_value=stat,widget.status
		endcase
 'QUIT'       : begin
                WIDGET_CONTROL,widget.base,/DESTROY
                endcase
  else: MESSAGE,"Event User Value Not Found.."
endcase
end

;---------------------------------------------------------------------------------------------------------------------------------------------
;Program that defines the GUI window and related structures....

pro full_tie
compile_opt idl2

close,/all

;common fd_data,fid
common widget_common, widget
common data_common, data

wsi=1024
;temperory dimension
;dim=unflip.dim
;fid=[unflip,flip]

;Widget Structure to hold all widget related information
widget = {widgetstruct, base:long(0), draw:long(0), status:long(0),rotang:long(0), loadseries:long(0), raw2align:long(0), sigma:long(0), $
          range:long(0),thrs:long(0), goalign:long(0), imnm:long(0), decon:long(0), mkmask:long(0), goDoGGT:long(0), goDoGLT:long(0), $
          shoimg:long(0), showimgnum:long(0), doseries:long(0), showsum:long(0), dotfs:long(0), revset:long(0), nummask:long(0),$
          showseries:long(0), shx:long(0), shy:long(0), dotierecon:long(0), nparams:long(0), method:long(0),showtfs:long(0),apptfs:long(0),$
          EE:long(0),CS:long(0),tiesetting:long(0),tiemethod:long(0),tikval:long(0),deriv:long(0),defstep:long(0),saveproj:long(0),loadproj:long(0),$
          dims:long(0),delta:long(0),numtfs:long(0),stype:long(0),wsi:long(wsi),drawid:long(0)}

;Data structure to hold common data processing parameters..
data = {datastruc, status:'Waiting for input..', revset:fix(0), rotang:fix(0), imnm:fix(0),flipflapfile:'flipunlfip.align',$
	params:fltarr(7),method:fix(2),click:intarr(25,2,2), curcorner:fix(1), gtlt:fix(1), nummask:fix(1), shx:fix(0), shy:fix(0),nparams:fix(4),$
	fl_exists:byte(0),unfl_exists:byte(0),dim:long(0),$
	ee:long(200000),cs:float(1.0e6),tiesetting:intarr(2),tikval:float(0.0),tiemethod:long(0),deriv:long(0)}


;Load the in-focus images from each set into the data structure
;openr,1,unflip.path+unflip.prefix+'.aligned'
;q=assoc(1,fltarr(dim,dim),12)
;data.unmask=q[unflip.num]
;data.un_img=q[0]*data.unmask
;close,1
;openr,1,flip.path+flip.prefix+'.aligned'
;q=assoc(1,fltarr(dim,dim),12)
;data.flmask=q[flip.num]
;data.fl_img=q[0]*data.flmask
;close,1
;
IF (XRegistered("full_tie") NE 0) then return

blcksz=350
DEVICE, GET_SCREEN_SIZE=scr
ydim = fix(scr[1] * 0.8) < wsi; set the viewing box to be 80% of the screen height
xdim = fix(ydim * 1.2) < wsi+blcksz ; make it a rectangular widget with aspect ratio 1.45
;ydim=wsi
;xdim=1450
;print,xdim,ydim

widget.base= WIDGET_BASE(TITLE='Full TIE Phase Retrieval', /COLUMN, SCR_XSIZE=xdim, SCR_YSIZE=ydim,x_scroll_size=xdim*0.9,y_scroll_size=ydim*0.9,xoffset=100,yoffset=300)
base2 = WIDGET_BASE(widget.base,/ROW)
base3 = WIDGET_BASE(base2,/column,xsize=blcksz)
base4 = WIDGET_BASE(base2,/column)
;info blocks

; create the drawing area.....
widget.draw= WIDGET_DRAW(base4,COLOR_MODEL=2,RETAIN=2,X_SCROLL_SIZE=ydim*0.8,Y_SCROLL_SIZE=ydim*0.8,/FRAME,/SCROLL,XSIZE=wsi, YSIZE=wsi)
widget.status= WIDGET_TEXT(base4,VALUE=string(" ",FORMAT='(A)'),XSIZE=20,UVALUE='STAT')

;action buttons and info start here....
infoblock = WIDGET_BASE(base3,/COLUMN,/FRAME,xsize=blcksz)

;Image series and loading/displaying image series
infoblock1= WIDGET_BASE(infoblock,/ROW)
descrip1 = WIDGET_BASE(infoblock1,/COLUMN,/BASE_ALIGN_LEFT)
fields1 = WIDGET_BASE(infoblock1,/COLUMN,/BASE_ALIGN_RIGHT)
; image series to consider;
item= WIDGET_LABEL(descrip1,VALUE='IMAGE SERIES',font=fontstrsmall,YSIZE=30,XSIZE=100,/ALIGN_LEFT)
widget.imnm = CW_BGROUP(fields1,['Unflip','Flip'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='IMNM',/NO_RELEASE)
infoblock1a=WIDGET_BASE(infoblock,/ROW,/ALIGN_LEFT)
;load series..
widget.loadseries = WIDGET_BUTTON(infoblock1a,/ALIGN_CENTER,VALUE='Load series',UVALUE='LOADSERIES',/FRAME,SENSITIVE=1)
widget.raw2align = WIDGET_BUTTON(infoblock1a,/ALIGN_CENTER,VALUE='Raw to Aligned',UVALUE='RAW2ALIGN',/FRAME,SENSITIVE=0)
widget.decon = WIDGET_BUTTON(infoblock1a,/ALIGN_CENTER,VALUE='Deconvolve with iMTF',UVALUE='DECONV',/FRAME,SENSITIVE=0)
;Image info..
infoblock1b = WIDGET_BASE(infoblock,/ROW)
descrip1b = WIDGET_BASE(infoblock1b,/COLUMN,/BASE_ALIGN_LEFT)
fields1b = WIDGET_BASE(infoblock1b,/COLUMN,/BASE_ALIGN_RIGHT)
lab = WIDGET_LABEL(descrip1b,VALUE='Dimensions',font=fontstrsmall,ysize=30)
widget.dims = WIDGET_TEXT(fields1b,VALUE=string(0,format='(I4)'),UVALUE='DIMS')
lab = WIDGET_LABEL(descrip1b,VALUE='Delta (nm/px)',font=fontstrsmall,ysize=30)
widget.delta = WIDGET_TEXT(fields1b,VALUE=string(1.0,format='(F10.2)'),UVALUE='DELTA',/EDITABLE)
lab = WIDGET_LABEL(descrip1b,VALUE='No. of TFS images',ysize=30,font=fontstrsmall)
widget.numtfs = WIDGET_TEXT(fields1b,VALUE=string(0,format='(I2)'),UVALUE='NUMTFS')
lab = WIDGET_LABEL(descrip1b,VALUE='Displayed Image #',ysize=30,font=fontstrsmall)
widget.showimgnum = WIDGET_TEXT(fields1b,/ALIGN_CENTER,VALUE='0',UVALUE='SHOWIMGNM',/EDITABLE)
infoblock1c = WIDGET_BASE(infoblock,/ROW)
seriestype=['Manual','Linear','Quadratic']
widget.stype=CW_BGROUP(infoblock1c,seriestype,UVALUE='SERIESTYPE',/exclusive,column=3,set_value=0,/NO_RELEASE)
;
;
;
;Image parameters for DoG setting/ ALignment for TFS and flip/unflip
infoblock2=WIDGET_BASE(base3,/COLUMN,/FRAME)
lab=WIDGET_LABEL(infoblock2,VALUE='ALIGNMENT PARAMETERS',/ALIGN_CENTER)
infoblock2b=WIDGET_BASE(infoblock2,/ROW)
descrip2b=WIDGET_BASE(infoblock2b,/COLUMN,/BASE_ALIGN_LEFT)
fields2b=WIDGET_BASE(infoblock2b,/COLUMN,/BASE_ALIGN_RIGHT)
; sigma
item= WIDGET_LABEL(descrip2b,VALUE='Sigma',font=fontstrsmall,YSIZE=30)
widget.sigma= WIDGET_TEXT(fields2b,VALUE=string(1,FORMAT='(F8.2)'),XSIZE=10,UVALUE='SIGMA',/EDITABLE)
; threshold parameter
item=WIDGET_LABEL(descrip2b,VALUE='Range',font=fontstrsmall,ysize=30)
widget.range= WIDGET_TEXT(fields2b,VALUE=string(0,FORMAT='(A)'),XSIZE=30,UVALUE='RANGE')
item=WIDGET_LABEL(descrip2b,VALUE='Threshold',font=fontstrsmall,ysize=30)
widget.thrs= WIDGET_TEXT(fields2b,VALUE=string(0,FORMAT='(F8.4)'),XSIZE=8,UVALUE='THRS',/EDITABLE)
; button to display the DoG filtered image
buttons2 = WIDGET_BASE(fields2b,/ALIGN_RIGHT,/ROW)
widget.goDoGGT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG > THR',UVALUE='DOGDISPLAYGT',/FRAME,SENSITIVE=1)
widget.goDoGLT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG < THR',UVALUE='DOGDISPLAYLT',/FRAME,SENSITIVE=1)
infoblock2c= WIDGET_BASE(infoblock2,/ROW)
descrip2c = WIDGET_BASE(infoblock2c,/COLUMN,/BASE_ALIGN_LEFT)
fields2c = WIDGET_BASE(infoblock2c,/COLUMN,/BASE_ALIGN_RIGHT)
item=widget_label(descrip2c,VALUE='No. of parameters',font=fontstrsmall,/ALIGN_LEFT,ysize=30)
;number of parameters for alignment
nparams = widget_base(fields2c,/align_center,/column)
widget.nparams = CW_BGROUP(nparams,['4','7'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='NPARAMS',/NO_RELEASE)
method = widget_base(fields2c,/align_center,/column)
;method of alignment
item=widget_label(descrip2c,VALUE='Method for Alignement',font=fontstrsmall,/ALIGN_LEFT,ysize=30)
widget.method = CW_BGROUP(method,['DoG','M.I'],/EXCLUSIVE,COLUMN=2,SET_VALUE=1,UVALUE='METHOD',/NO_RELEASE)
infoblock2d=WIDGET_BASE(infoblock2,/ROW)
descrip2d=WIDGET_BASE(infoblock2d,/ROW,/ALIGN_LEFT)
;Mask for alignment
item= WIDGET_LABEL(descrip2d,VALUE='Define Rectangular Gaussian Mask',font=fontstrsmall,/ALIGN_CENTER)
widget.nummask= WIDGET_TEXT(descrip2d,VALUE=string(0,FORMAT='(I3)'),XSIZE=5,UVALUE='NUMMASK',/EDITABLE,/ALIGN_CENTER)
widget.mkmask = WIDGET_BUTTON(descrip2d,/ALIGN_CENTER,VALUE='Go',UVALUE='MKMASK',/FRAME,SENSITIVE=1)
;
;
;
;TFS alignment 
infoblock3=WIDGET_BASE(base3,/COLUMN,/FRAME)
lab=WIDGET_LABEL(infoblock3,VALUE='TFS ALIGNMENT',/ALIGN_CENTER)
butts3=WIDGET_BASE(infoblock3,/ROW,/ALIGN_CENTER)
widget.dotfs=WIDGET_BUTTON(butts3,VALUE='Do Align',UVALUE='DOTFSALIGN',/ALIGN_CENTER,/FRAME,SENSITIVE=0)
widget.showtfs=WIDGET_BUTTON(butts3,VALUE='Show Align',UVALUE='SHOWTFSALIGN',/ALIGN_CENTER,/FRAME,SENSITIVE=0)
widget.apptfs=WIDGET_BUTTON(butts3,VALUE='Apply Align',UVALUE='GOTFSALIGN',/ALIGN_CENTER,/FRAME,SENSITIVE=0)
;
;
;
;Parameter setting for unflip/flip align
infoblock4=WIDGET_BASE(base3,/COLUMN,/FRAME)
lab=WIDGET_LABEL(infoblock4,VALUE='UNLIP/FLIP ALIGN',/ALIGN_CENTER)
infoblock4a= WIDGET_BASE(infoblock4,/ROW)
descrip4a = WIDGET_BASE(infoblock4a,/COLUMN,/BASE_ALIGN_LEFT)
fields4a = WIDGET_BASE(infoblock4a,/COLUMN,/BASE_ALIGN_RIGHT)
; overall shift parameters to correct for a poorly positioned image
item= WIDGET_LABEL(descrip4a,VALUE='Image Shift',font=fontstrsmall,YSIZE=30)
buttons3 = WIDGET_BASE(fields4a,/ALIGN_CENTER,/ROW)
widget.shx= WIDGET_TEXT(buttons3,VALUE=string(0,FORMAT='(I5)'),XSIZE=8,UVALUE='SHX',/EDITABLE)
widget.shy= WIDGET_TEXT(buttons3,VALUE=string(0,FORMAT='(I5)'),XSIZE=8,UVALUE='SHY',/EDITABLE)
item= WIDGET_LABEL(descrip4a,VALUE='Rotation Angle',font=fontstrsmall,YSIZE=30,/ALIGN_CENTER)
;buttons3 = WIDGET_BASE(fields3b,/ALIGN_CENTER,/ROW)
; overall rotate and reverse setting parameters
widget.rotang= WIDGET_TEXT(fields4a,VALUE=string(0,FORMAT='(I5)'),XSIZE=8,UVALUE='ROTANG',/EDITABLE)
item = WIDGET_LABEL(descrip4a,VALUE='REVERSE SET',font=fontstrsmall,YSIZE=30)
rev_set = WIDGET_BASE(fields4a,/align_right,/column)
widget.revset= CW_BGROUP(rev_set,['0','1'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='REVSET',/NO_RELEASE)
butts = WIDGET_BASE(infoblock4,/align_center,/row)
widget.showsum = WIDGET_BUTTON(butts,/ALIGN_CENTER,VALUE='Show Sum',UVALUE='SHOWSUM',/FRAME)
; buttons to start the alingment and loop through after alingment
widget.doseries= WIDGET_BUTTON(butts,/ALIGN_CENTER,VALUE='Do Align',UVALUE='DOALIGN',/FRAME,SENSITIVE=0)
widget.showseries = WIDGET_BUTTON(butts,/ALIGN_CENTER,VALUE='Show Align',UVALUE='SHOWALIGN',/FRAME,SENSITIVE=0)
widget.goalign= WIDGET_BUTTON(butts,/ALIGN_CENTER,VALUE='Apply Align',UVALUE='GOALIGN',/FRAME,SENSITIVE=0)
;
;
;
;Parameters for TIE Reconstruction....
infoblock5= WIDGET_BASE(base3,/COLUMN,/FRAME)
lab=WIDGET_LABEL(infoblock5,VALUE='TIE RECON',/ALIGN_CENTER)
infoblock5b = WIDGET_BASE(infoblock5,/ROW)
descrip5b = WIDGET_BASE(infoblock5b,/COLUMN,/BASE_ALIGN_LEFT)
fields5b = WIDGET_BASE(infoblock5b,/COLUMN,/BASE_ALIGN_RIGHT)
labelEE=WIDGET_LABEL(descrip5b,VALUE='EE (V) :',ysize=30)
widget.EE=WIDGET_TEXT(fields5b,VALUE=string(200000,format='(F10.2)'),$
	UVALUE='EE',/EDITABLE)
labelcs=WIDGET_LABEL(descrip5b,VALUE='Cs (nm) :',ysize=30)
widget.cs=WIDGET_TEXT(fields5b,VALUE=string(1.0e6,format='(E10.2)'),$
        UVALUE='CS',/EDITABLE)
labeldefstep=WIDGET_LABEL(descrip5b,VALUE='Defocus step (nm) :')
widget.defstep=WIDGET_TEXT(fields5b,VALUE=string(0,format='(F10.2)'),$
        UVALUE='DEFSTEP')
infoblock5c = WIDGET_BASE(infoblock5,/ROW)
settingvals=['Tikhonov','Symmetrize']
widget.tiesetting=CW_BGROUP(infoblock5c,settingvals,UVALUE='TIESETTING',/nonexclusive,column=2,set_value=[1,1])
widget.tikval=WIDGET_TEXT(infoblock5c,VALUE=string(2.0,format='(F6.2)'),UVALUE='TIKVAL',/EDITABLE)
methods=['Laplacian','Inverse Gradient']
widget.tiemethod=CW_BGROUP(infoblock5,methods,UVALUE='METHODS',/exclusive,column=2,set_value=0)
derivmethods=['3-point dI/dz','PolyFit dI/dz']
widget.deriv=CW_BGROUP(infoblock5,derivmethods,UVALUE='DERIVMETHOD',/exclusive,column=2,set_value=0,/NO_RELEASE)
widget.dotierecon=WIDGET_BUTTON(infoblock5,/ALIGN_CENTER,VALUE='Do TIErecon',UVALUE='DOTIERECON',/FRAME,SENSITIVE=0)
;
;End
buttons5 = WIDGET_BASE(base3,/ALIGN_CENTER,/ROW)
widget.saveproj = WIDGET_BUTTON(buttons5,/ALIGN_CENTER,VALUE='Save Project',UVALUE='SAVEPROJ',/FRAME,SENSITIVE=0)
widget.loadproj = WIDGET_BUTTON(buttons5,/ALIGN_CENTER,VALUE='Load Project',UVALUE='LOADPROJ',/FRAME,SENSITIVE=0)
quit  = WIDGET_BUTTON(buttons5,/ALIGN_CENTER,VALUE='Quit',UVALUE='QUIT',/FRAME)
;
WIDGET_CONTROL, /REALIZE, widget.base
WIDGET_CONTROL, widget.draw, GET_VALUE=drawID
widget.drawID = drawID
;
;check if alignment file exists and update numbers accordingly
;fi=file_info(data.alignfile)
;if (fi.exists eq 1) then begin
;  openr,1,data.alignfile
;  tmp=0
;  readf,1,tmp
;  if (tmp eq 4) then begin
;    params=fltarr(4)
;    WIDGET_CONTROL,SET_VALUE=0,widget.nparams
;    data.nparams=4
;  endif else begin
;    params=fltarr(7)
;    WIDGET_CONTROL,SET_VALUE=1,widget.nparams
;    data.nparams=7
;  endelse
;  al_set=fltarr(4)
;  readf,1,al_set
;  readf,1,params
;  close,1
;  if (data.nparams eq 4) then data.params[0:3]=params else data.params[0:6]=params
;  data.shx=al_set[0]
;  data.shy=al_set[1]
;  data.rotang=al_set[2]
;  data.revset=al_set[3]
;  WIDGET_CONTROL,SET_VALUE=string(al_set[0],format='(I5)'),widget.shx
;  WIDGET_CONTROL,SET_VALUE=string(al_set[1],format='(I5)'),widget.shy
;  WIDGET_CONTROL,SET_VALUE=string(al_set[2],format='(I5)'),widget.rotang
;  WIDGET_CONTROL,SET_VALUE=fix(al_set[3]),widget.revset
;  WIDGET_CONTROL,SENSITIVE=1,widget.showseries
;  WIDGET_CONTROL,SENSITIVE=1,widget.goalign
;  WIDGET_CONTROL,SET_VALUE='Alignment file exists. Pre-loaded the values. Press Show Align to check.',widget.status
  ;Apply the alignment parameters
;  if (data.revset eq 1) then begin
;	 imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
;  endif else begin
;         imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
;  endelse
;  imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
;  line = findgen(data.dim)-float(data.dim/2)
;  bintx = replicate(1.0,data.dim)##line
;  binty = replicate(1.0,data.dim)#line
;  if (data.nparams eq 4) then begin
;      imb = rot(imb,data.params[2],data.params[3],cubic=-0.5,missing=0.0)
;      imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+params[1],cubic=-0.5,missing=0.0)
;  endif else begin
;      imb = rot(imb,data.params[2],cubic=-0.5,missing=0.0)
;      imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
;			bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=0.0,cubic=-0.5)
;  endelse
;  data.fl_img_al = imb
;endif
;
; Make the draw widget the current IDL drawable area.
WSET, drawID
;tvscl,rebin(data.un_img,widget.wsi,widget.wsi)
;
XMANAGER,"full_tie",widget.base,EVENT_HANDLER="full_tie_ev",/NO_BLOCK
;
;
return
end

