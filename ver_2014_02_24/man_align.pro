function man_align,imstack=imstack,params=params,verbose=verbose
compile_opt idl2
;+
;program to manually a set of 3 images from TFS.
;assumes the imstack=fltarr(dim,dim,3) and aligns
;0 & 2 wrt 1. Returns the params=[shiftx,shifty,rot,mag]
;@author, Charudatta Phatak, CMU, 02/13/09.
;-

sz=size(imstack,/dim)
dim=sz[0]
num=sz[2]
params=fltarr(4,num)
params[0:3,1]=[0.0,0.0,0.0,1.0]
wsi=1024

;define keyboard sequences
enter=total(byte(10))
up=total(byte([27,91,65]))
down=total(byte([27,91,66]))
left=total(byte([27,91,68]))
right=total(byte([27,91,67]))
shup=total(byte([27,91,49,59,50,65]))
shdown=total(byte([27,91,49,59,50,66]))
shleft=total(byte([27,91,49,59,50,68]))
shright=total(byte([27,91,49,59,50,67]))
plus=total(byte(45))
minus=total(byte(43))
space=total(byte(32))
quit=total(byte(113))
rotr=total(byte(114))
rotl=total(byte(108))

;Print out the help message to tell how to modify the images..
Print,'-----------Manual Alignment Help---------------'
print,'Arrow keys move the image by 1 pixel'
print,'Shift+arrow move the image by 0.5 pixel'
print,'+/- changes the magnification of the image'
print,'R/L rotates the image to Right/Left by 1 deg'
print,'space - flicks the two images 10 times.'
print,'Press enter when done aligning the current images.'
print,'q - quits the manual alignment at any time.'
print,'-----------------------------------------------'

im1 = reform(imstack[*,*,0])
window,5,xsi=wsi,ysi=wsi,retain=2,title='Manual alignment window'
for nim=0,num-1,2 do begin
	if (nim eq 0) then begin
		print,'Aligning Underfocus image wrt to Infocus image..'
		im2 = reform(imstack[*,*,1])
	endif else begin
		print,'Aligning OVerFocus image wrt to Infocus images..'
		im2 = reform(imstack[*,*,2])
	endelse
	flag=0
	nim2=im2
	tvscl,rebin(0.5*(im1+nim2),wsi,wsi)
	bintx=findgen(dim)#replicate(1.0,dim)
	binty=replicate(1.0,dim)#findgen(dim)
	xs=0.0 & ys=0.0 & rot_ang=0.0 & mag=1.0
	while(flag eq 0) do begin
		wset,5
		key=total(byte(get_kbrd(/escape)))
		case key of
			enter : begin
				flag=1
				end
			up : begin
			     binty-=1.0
			     ys-=1.0
		     	     end
			down: begin
			      binty+=1.0
			      ys+=1.0
		      	      end
			left: begin
			      bintx+=1.0
			      xs+=1.0
		      	      end
			right: begin
			       bintx-=1.0
			       xs-=1.0
			       end
		        shup: begin
			      binty-=0.5
			      ys-=0.5
			      end
			shdown: begin
			        binty+=0.5
			        ys+=0.5
			        end
			shleft: begin
			        bintx+=0.5
			        xs+=0.5
			        end
			shright: begin
			         bintx-=0.5
			         xs-=0.5
			         end
			plus: begin
				mag +=0.01
				end
			minus: begin
				mag -=0.01
				end
			rotr: begin
				rot_ang += 1.0
				end
			rotl: begin
				rot_ang -=1.0
				end
			space: begin
				for icnt=0,10 do begin
					tvscl,rebin(im1,wsi,wsi)
					wait,0.3
					tvscl,rebin(nim2,wsi,wsi)
					wait,0.3
				endfor
				end
			quit: begin
				wdelete,5
				return,-1
			      end
		else : print,'Invalid key'
		endcase
		nim2=rot(im2,rot_ang,mag,cubic=-0.5,missing=0.0)
		nim2=interpolate(nim2,bintx,binty,cubic=-0.5,missing=0.0)
		tvscl,rebin(0.5*(im1+nim2),wsi,wsi)
		xyouts,0.05,0.05,string(total((im1-nim2)^2)),/normal,color=200,charsize=1.5
	endwhile
	params[0:3,nim] = [xs,ys,rot_ang,mag]
	print,'Final Alignment Params:',params[0:3,nim]
endfor
wdelete,5
return,1
end

