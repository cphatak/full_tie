@read_dm3
pro rd_fls,file_data
;
;this routine is similar to rd_fls in through focus series folder
;written by MDG. Essentially this one only reads in all the image files
;and returns a str. containing the information.
;
;line 1: n = no. of images
;line 2 - n+1 : filenames
;-
path="/Users/cphatak/ANL_work/"
;path='/Volumes/viz/Labs/ALTEM/CD/2013_Dec_11/tomo_CFB_Py_SiO'
;path="/Volumes/data/CD_postdoc/TEM/ALTEM/2012/020812/"
fname=dialog_pickfile(title='Select *.fls file',/noconfirm,path=path,$
	get_path=fpth,filter='*.fls')
if (fname eq '') then begin $
		file_data=-1
	goto, cancel
endif

;read the fls file
rlen=strlen(fname)
plen=strlen(fpth)
prefix=strmid(fname,plen,rlen-plen-4)
print,'file name prefix = ',prefix
openr,1,fname
num=0
dim=0
readf,1,num
ds={num:num, $
	prefix:prefix, $
	names:strarr(num), $
	path:fpth, $
	fexists:bytarr(num), $
	dim:dim,$
	delta:0.0}
print,'No. of images detected ',num
nm=''
for i=0,num-1 do begin
	readf,1,nm
	ds.names[i]=nm
	fi=file_info(fpth+nm)
	if (fi.exists eq 1) then begin
		ds.fexists[i]=1B
		print,'Checking for file ',fpth+nm,'-> OK'
	endif else begin
		ds.fexists[i]=0B
		print,'Checking for file',fpth+nm,'does not exist.'
	endelse
endfor
close,1
if (total(ds.fexists) ne byte(num)) then begin
	print,'One or more files are missing. Please check.'
	goto, cancel
endif

;getting dimensions
read_dm3,ds.path+ds.names[0],iminfo;,/verbose
imstack=iminfo.image
delta=iminfo.scale
sz=iminfo.dims

;dims=size(image,/dim)
print,'Image dimensions are ',sz
ds.dim = sz[0]
ds.delta = delta[0]

problem=0
if (total(ds.fexists) ne ds.num) then begin
	print,'Error: One of more files not found.'
	problem=1
endif

if problem eq 0 then file_data=ds else file_data=-1
cancel:
end
