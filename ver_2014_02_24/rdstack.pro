@read_dm3
;+
;NAME: RDSTACK.pro
;
;USAGE: rdstack, fs
;
;PURPOSE: The routine reads in a stack of dm3 images, recorded using the
;         AcqFocalSeries_CD.s script in DM, and returns a pointer with the
;         information regarding the stack.
;
;ARGUMENT: fs - pointer to hold the information structure for the file
;
;REQUIRES: READ_DM3.pro - to read the dm3 file directly.
;
;AUTHOR: CD Phatak, ANL, 09/14/2012.
;-
pro rdstack,fs
compile_opt idl2
;Set the initial path to look at - 
path="/Users/cphatak/ANL_work/";VFET/"
;path="/Volumes/DATA2_RAID/CD_work/TEM/"

;Pick the stack file.
fname=dialog_pickfile(title='Select dm3 stack file',/noconfirm,path=path,$
	get_path=fpth,filter='*.dm3')
if (fname eq '') then begin $
		file_data=-1
	goto, cancel
endif

;get filename,pathname from the file
rlen=strlen(fname)
plen=strlen(fpth)
prefix=strmid(fname,plen,rlen-plen-4)
print,'file name prefix = ',prefix

;create file structure for storing data
ds={num:0, $
	prefix:prefix, $
	path:fpth, $
	dim:0L, $
  delta:1.0, $
  defstep:0.0}

;read the stack and get the info
read_dm3,ds.path+ds.prefix+'.dm3',iminfo;,/verbose
imstack=iminfo.image
delta=iminfo.scale
sz=iminfo.dims
sz=size(imstack,/dim)
print,'Stack Dimensions:',sz
if (n_elements(sz) ne 3) then begin
  print,'Image Stack not 3D.'
  goto, cancel
endif
dimx = sz[0] & dimy = sz[1]
if (dimx ne dimy) then begin
  print,'Stack does not contain square images.'
  goto, cancel
endif
ds.dim = dimx
ds.num = sz[2]
ds.delta = delta[0]
ds.defstep = delta[2]
;return the file structure
fs = ds

cancel:
end
