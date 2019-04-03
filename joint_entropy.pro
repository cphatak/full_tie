function joint_entropy,images,mutual_information=mutual_information,feature_space=feature_space,NMI=NMI,ECC=ECC,help=help,nozero=nozero
;
; this function returns the joint entropy, and optionally the feature space and/or
; the mutual information of two grayscale images 
;
; arguments:
; images               : the image bytarr[num_of_images,cols,rows] (must be gray scale)
;
; options:
; /help lists the arguments and options
; mutual_information   : stores the conventional mutual information in a named variable
; /NMI                 : stores the normalized mutual information instead
; /ECC                 : stores the Entropy Correlation Coefficient instead
; feature_space        : returns the 2D feature space (i.e. joint histogram) in a named variable
;
if keyword_set(help) then begin
	print,'Here are the arguments and options for the joint_entropy routine'
	print,'arguments:'
	print,'	images               : the image bytarr[num_of_images,cols,rows] (must be gray scale)'
	print,' '
	print,'options:'
	print,'	/help lists the arguments and options and exits without computation'
	print,'	mutual_information   : stores the conventional mutual information in a named variable'
	print,'	/NMI                 : stores the normalized mutual information instead'
	print,'	/ECC                 : stores the Entropy Correlation Coefficient instead'
	print,'	feature_space        : returns the 2D feature space in a named variable'
	print,' '
	return,0
endif

; get the feature space for these two images
fs = float(joint_histogram(images))
if keyword_set(nozero) then fs[0,0] = 0.0

; normalize feature space by the total number of entries
t=total(fs)

; compute the joint entropy
q=where(fs gt 0.0,cnt)
je_AB = -total(fs[q]/t*alog(fs[q]/t))

; export feature space if requested
if arg_present(feature_space) then feature_space = fs

; compute mutual information in one of three forms
if arg_present(mutual_information) then begin
	; entropy of image A
	hA = total(fs/t,1)       ; sum the feature space onto the first axis
	hA = hA/total(hA)
	q = where(hA gt 0.0,cnt)
	hA = -total(hA[q]*alog(hA[q]))
	; entropy of image B
	hB = total(fs/t,2)       ; sum the feature space onto the second axis
	hB = hB/total(hB)
	q = where(hB gt 0.0,cnt)
	hB = -total(hB[q]*alog(hB[q]))
	; standard mutual information
	mutual_information = hA + hB - je_AB
	; normalized mutual information
	if keyword_set(NMI) and (je_AB ne 0.0) then mutual_information = (hA+hB)/je_AB
	; entropy correlation coefficient
	if keyword_set(ECC) then mutual_information = 2.0 - 2.0*je_AB/(hA+hB)
endif

; that's it...
return,je_AB
end