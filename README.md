full_tie
========

This repository contains codes for performing Full Transport-of-Intensity Equation based phase retrieval using IDL

Wiki for FULL_TIE.PRO

C. Phatak
April 3, 2019.

Running Environment:
This repository contains all the codes for performing full TIE based phase reconstruction. It requires two through-focus series of images: one with the sample as-is (unflip) and one with the sample flipped upside down (flip). Theses codes run in IDL environment. If you have a IDL license that is great. If not, you can download the IDL virtual machine which is available for free (https://www.harrisgeospatial.com/Software-Technology/IDL). Then run the IDL virtual machine and load the full_tie.sav file which should load all the functions necessary to run the program. If you have any issues please contact me at c.phatak@gmail.com

Instructions for using the program:
1) The data must be input as either an fls file or a DM3 3D image stack. If it is a 3D image stack in DM3 format, the program will automatically read the defocus values from the DM3 header as well as the pixel resolution. 
The *.fls file is a simple text file and can be created using any texteditor. The format of the file should be as follows:


    Number of images in series  
    Filenames of all the images starting with infocus.dm3  
    -defocus.dm3 (underfocus in increasing defocus value)  
    …  
    +defocus.dm3 (overfocus in increasing defocus values)  
    defstep  
    2 * defstep  
    ....  
    n/2 * defstep (highest defocus)  


Save the fls file with a useful name since that will be the prefix that will used to save all the subsequent images and intermediate files. Select the appropriate series before loading at the top – “unflip” or “flip”. The rest of the settings are applied accordingly to that series. 

2) After the images are loaded, you should be able to see the information about them populated in the left hand side such as dimensions, no. of tfs images etc. Then select the type of through focus series such as Manual, Linear or Quadratic series from the radio buttons. This will help compute the defocus values for respective images.

3) The program can then be used for aligning the through focus series images. Set the alignment parameters such as the sigma, threshold, no. of parameters, method of alignment and # of mask for alignment. If the no of parameters selected is 4, the alignment is done to account for uniform rotation and scaling, it the no. of parameters is 7, then anisotropic scaling(affine) is also accounted for. The method of alignment can be chosen from Mutual information based (MI) or difference of gaussian (DOG). For DoG method, binary versions of the images are used for aligning. The binary images are produced by changing the value of sigma (the sigma for gaussians) and threshold value. Everytime the sigma is changed, the intensity range is updated so that the threshold value can be selected appropriately. A lot of times this is more accurate but very slow. The MI method uses mutual information between the two images to align them and is quite fast. Typically for ‘4’ parameter alignment and for ‘7’ parameter with very low affine corrections, the MI method will work the best and quickest. Once the alignment is complete, click show align to see the aligned images. If the alignment is good, then click apply align so that the alignment is actually applied to the image stack and a *.aligned file will be created. 

4) If the stack alignment was done using using the ImageJ plugin “stackreg”, then the output file should be saved as a raw file with the same prefix as the fls file with the added extension at the end ‘_alimj.raw’. Then click the ‘Raw to Aligned’ button. It should then load the aligned stack and go through the stack back and forth in the image window. And save the aligned image stack in a *.aligned file. 

5) Once the *.aligned file is created, click the “Deconvolve with iMTF” button to deconvolve the the modulation transfer function from the images. The deconvolution will occur automatically and all the image stack will be saved as a *.decon file. The *.decon file will be the one used for further processing such as alignment of flip/unflip series and finally the phase reconstruction.

6) Repeat these steps with the other series. So that now there should be *.decon files for both the “unflip” and “flip” series.

7) Once that is complete, the unflip and flip data sets need to be aligned wrt each other. Click the Show Sum button to show the sum of unfliped infocus image and flipped infocus images. The settings that can be updated initially are the reverse_set (to reverse the image flip), rotation angle (to correct the rotation of the image) and image shift (manual image shifts in x & y). The goal is using these initial settings to get the images to lie on top of each other as close as possible. The final subpixel alignment will be done automatically. Once the manual values are good, click on the Do Align button to perform the automated alignment. The alignment parameters that will be used will be the same from the “Alignment Parameters” window. They can be changed for this specific alignment as needed. In this case, typically the DoG method gives the best alignment results although quite slow. 

8) Once the alignment is complete, the alignment parameters will be saved in a file. Click on “Show align” in the UNFLIP/FLIP ALIGN box to see the alignment. This will flicker the two images back and forth. If you are satisfied with the alignment, then click on “Apply Align” to apply the alignment to the flipped stack of images. 

9) After all the alignment is complete, now we can perform the phase reconstruction. The first three boxes in the “TIE RECON” box have pre-populated values for accelerating voltage, E (in V), Cs (in nm) and defocus. You can change the E and Cs as per the microscope. The defocus will automatically be used as per the previous definition using manual, linear or quadratic series. Next setting is for symmetrization  of images and using Tikhonov filtering. Those can be set “ON” or “OFF” by pressing the buttons next to them. The box with the numerical value is the value for the tikhonov filtering. NOTE: FOR NON-ZERO VALUES OF TIKHONOV FILTER, THE PHASE RECONSTRUCTION IS NOT QUANTITATIVE.
The next set of options is the method for TIE reconstruction. You can either select “Laplacian” method or the “Inverse gradient” method for solving the TIE. The final option is for computing the intensity derivative. You can either choose the central difference method for intensity derivative or a polynomial fit based intensity derivative calculation. After all the settings are done, click on “Do TIErecon” to start the phase reconstruction. The phase reconstruction uses the updated method of separately reconstructing the magnetic and electrostatic phase shifts using a modified definition for intensity derivative.vi The program will ask whether to resize the entire image to 1024 by 1024 size or to select a 1K by 1K  region from the image. After the reconstruction is complete, a folder called “images” will be created in the same location as the flipped image stack, and images of reconstructed, Bx.t, By.t, Bb.t (magnitude), colormap, intensity derivatives, phase maps and infocus images will be saved. All the images (except colormap) are saved as 32-bit float files, which can be opened using ImageJ or DM for further analysis. The pixel calibration in the reconstructed images is not yet implemented correctly. 

References:
 (i) E. Gulsoy, J. Simmons, and M. De Graef, Scr. Mater. 60, 381 (2009).  
 (ii) http://bigwww.epfl.ch/thevenaz/stackreg/  
 (iii) V. V Volkov, Y. Zhu, and M. De Graef, Micron 33, 411 (2002).  
 (iv) K. Ishizuka and B. Allman, J. Electron Microsc. (Tokyo). 54, 191 (2005).  
 (v) M. De Graef, Introduction to Conventional Tranmission Electron Microscopy, Cambridge Press, 2003.  
 (vi) E. Humphrey, C. Phatak, A. K. Petford-Long, and M. De Graef, Ultramicroscopy 139, 5 (2014).  
