//An UI for a collection of drift alignment methods
//
//Useage:
//(Optional: place a rectangular ROI around a feature in the displayed slice of the image stack)
//Load 3D image stack via the "Assign" button
//
//"Hanning" places a Hanning window in the real-space domain in order to reduce FFT artefacts from the sharp image edges 
//	Creates a new copy in memory which can quickly overwhelm the PC, use accordingly
//
//Load existing aperture (image dimensions must match), or create new:
//"Top-hat" places a sharp aperture around the image center
//"B'worth" places a smooth aperture around the image center; order (slope) can be adjusted
//"Positive" places a positive mask around the brightest spot within the non-volatile ROI in Fourier Transform preview
//	Positive spot masks are invisible in the preview until the mask is saved or applied
//"Negative" suppresses the brightest spot in the FT ROI by multiplying a negative peak function
//	Useful in combination with the Buttersworth filter to eliminate residual sideband signal
//"Blur" applies a Gaussian convolution kernel to the mask
//
//"Sideband" performs an iterative maximum search around the brightest spot outside the center pixel, and attempts to remove them
//"Fresnel" places a Gauss-line filter over the Fresnel streak
//	Uses the Holoworks (tested on v5.12) function; disable or remove if this leads to conflicts
//
//All mask-forming steps can be reverted once via the "Undo" button
//
//"XCF" and "PCF" measure the drift via cross-correlation / phase-correlation methods and create a series of shift vectors dy/dx as output
//	For PCF references see:
//	C.D. Kuglin, D.C. Hines, Proceedings of the IEEE International Conference on Cybernetics and Society, 1975, p. 163
//	Meyer, R., Kirkland, A., Saxton, W., 2002. Ultramicroscopy 92, 89
//	T. Niermann, M. Lehmann, 2014, Micron 63, p. 28-34
//
//The "Pairwise" switch determines whether subsequent images or drift vs. first in the time series are considered
//"Max X" / "Max Y" fields limit the search range, e.g. to within less than once unit cell
//	Default search window size is 1/2 of the total image, but this is much larger than realistically applicable since the site of interest would very quickly move outside the camera frame
//
//"Load" (Drift) allows importing an existing drift measurement
//"Apply" aligns the slices according to the currently measured / loaded drift series
//
//By default, the field of view is cut to the image region that overlaps between all slices. Set use_cut = 0 to perform this manually later on.
//
//***All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.
//***temscript is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENCE.txt file for any details.
//***All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.

string version = "v1.0 Oct. 2023 U.Ross, J.Lindner"
//flags for debugging
number verbose = 1	//extra feedback in the results window
number use_cut = 1	//cut edges of the drift-corrected stack
number has_HW = 1	//access to Holoworks scripting functionality; if the script refuses to run without the Holoworks library, delete the contents of the void remove_fresnel(object self){ function

//the interface lets us use block/reset from anywhere within the script, not strictly necessary
interface ThreadInterface
{
	void block_UI(object self);
	void reset_UI(object self);
}

//the UI class, all functions are embedded
class DialogUI:uiframe{
//global vars
image temp, front, fft, oldfft, cfft, mask, oldmask	//persistent images
image posmask 										//positive mask image overrides other masks
number t,l,b,r										//ROI vertices
number sx,sy,sz,rx,ry								//dimensions of the input image stack and ROI
number hat_radius,maxx,maxy,pairs					//parameters to be read from the UI
number outID										//passing output between functions by reference
string log											//log string to retrace your steps, is added to the output image as description text

//various utility functions
//return ROI vertices left-top and right-bottom (note the order left-top-right-bottom instead of top-left-bottom-right due to personal preference)
number get_ROI_position(object self, image in, number &px, number &py, number &bx, number &by){
	imagedisplay imgdisp = imagegetimagedisplay(in,0)
	
	try	
	{
		ROI theROI = ImageDisplayGetROI( imgDisp, 0 )	//assuming only one ROI present
		ROIisvalid(theROI)
		ROIIsRectangle(theROI)
		ROIGetVertex( theROI, 0, px, py )
		ROIGetVertex( theROI, 2, bx, by )
	}
	catch
	{
		getsize(in,bx,by)
		px = 0; py = 0;	//return the original image vertices if no ROI found 
		if(verbose) result("\nROI invalid or does not exist, returning whole image");
		return(0)
	}
	return(1)	
}

//cross-wise parabolic fit, originally by D.G. Mitchell (dmscripting.com) and adjusted for our purposes
//returns by reference the absolute coordinates (relative to top-left origin) within the search window, i.e. the origin may need to be adjust later on
number ImageRefineExtrema(object self, image img, number &px, number &py){
		number tx1, tx2, ty1, ty2
		number numa, numb, numc
		number valX, valY

		// do a parabolic fit in direction x
		numa = img.GetPixel(px-1, py)
		numb = img.GetPixel(px  , py)
		numc = img.GetPixel(px+1, py)
		if(!(max(numa,numb) == numb) && (max(numb,numc) == numb)){ result("\nNot a local maximum."); return 0; }	//restrict to local maxima in x and y direction
		
		tx2=(numa-numc)/2.0
    	tx1=(numa+numc)/2.0-numb
		valX=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb

		// do a parabolic fit in direction y
		numa = img.GetPixel(px, py-1)
		numb = img.GetPixel(px, py  )
		numc = img.GetPixel(px, py+1)
		if(!(max(numa,numb) == numb) && (max(numb,numc) == numb)) return 0;
		
		ty2=(numa-numc)/2.0
    	ty1=(numa+numc)/2.0-numb
		valY=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb
		// update pixel position if the shift is smaller than one pixel; otherwise the coarse maximum search has failed / landed at the edge of the restricted search window
		number dx = tx2/(2.0*tx1)
		number dy = ty2/(2.0*ty1)
		if(( abs(dx) <= 1 ) && ( abs(dy) <= 1 )){
			px += dx;
			py += dy;
		}
		else{
			result("\nFit refinement outside 1 pixel deviation")
			return 0
		}
		
		// check to see whether to return mimimum or maximum; no longer applicable but left in in case we want to revert the functionality
		If(numa>=numb && numc>=numb) return min(min(valX, valY), numb)
		If(numa<=numb && numc<=numb) return max(max(valX, valY), numb)

		Return numb
}

//take an input mask and dampens it around a pixel relative to the image center
//multiplies a negative peak function onto the input image, centered around the coordinate of choice
//peak shape could be anything but I chose a flexible Pearson function,  in partice width and order are hard-coded
image Pearson(object self, image in, number std, number order, number jx, number jy){
	in = in * (1-( (std**(2*order))/( std**2+(2**(1/order)-1)*(sqrt((icol-jx)**2+(irow-jy)**2))**2)**order ))
	return in
}

//thread init
object init(object self, number threadID){
		//passing in the scriptobject ID which could be used to do sanity checking, in practice omitted
		if(verbose) result("\n\nUI interface to align image stacks. If you want to use a ROI, place it before loading the stack.");
		return self
	}

//button reponse functions
//load the front image stack and refresh all auxiliary variables to prevent dimension errors
void load_stack(object self){
		try{ deleteimage(temp); deleteimage(mask); deleteimage(oldmask); deleteimage(fft); deleteimage(oldfft); deleteimage(cfft); }
		catch break;
		
		front := getfrontimage()
		try front.get3dsize(sx,sy,sz);
		catch{ okdialog("Image is not a 3D stack!"); exit(0); }
		
		t=0; l=0; b=sy; r=sx; rx=sx; ry=sy;		//critical to fill in rx/ry
		log = ""
		
		number hasROI = self.get_ROI_position(front,l,t,r,b)	//test for ROI and read coordinates if present
		if(hasROI){
				rx = (r-l)
				ry = (b-t)
				temp = slice2(front,l,t,0,0,rx,1,1,ry,1)		//grabbing the ROI in the first slice, does not have to be square (!)
				cfft = realFFT(temp)							//cfft is complex-valued, which is why we treat it separately
				maxx = rx/4										//default search window guess for ROI is +/- 1/4th
				maxy = ry/4
				log += "\nLoaded stack "+getname(front)+" with ROI at "+l+"/"+t+"/"+r+"/"+b+""
			}
		else{
				cfft = realFFT(slice2(front,0,0,0,0,sx,1,1,sy,1))		//also only displaying the FFT of the first slice.
				maxx = sx/4												//default search window guess for full frame
				maxy = sy/4
				log += "\nLoaded stack "+getname(front)+" without ROI"
			}
		
		fft = modulus(cfft)
		fft = (fft-min(fft))/(max(fft)-min(fft))	//renormalization for the preview frame
		fft.showimage()
		hat_radius = min(rx/4,ry/4)					//default radius guess; the height (y-size) of the ellipse in case of non-square image/ROI
		DLGvalue(self.lookupelement("hatradius"),hat_radius)	//write the values to the dialog UI
		DLGvalue(self.lookupelement("maxx"),maxx)
		DLGvalue(self.lookupelement("maxy"),maxy)
		mask := realimage("",4,rx,ry)+1				//default mask is transparent (==1)
		mask.setname("FFT aperture")
		oldmask = mask								//oldmask retains the previous mask for the undo button
		oldfft = fft								//oldfft contains a permanent copy of modulus(FT)
		posmask = mask*0							//positive masks are not fully implemented, by default opaque (==0)
		
		//create a nonvolatile ROI in the FFT display to locate maxima; size is abitrary since only the center will matter
		ROI theROI = NewROI()
		ROISetRectangle( theROI, (ry/2)-(ry/128), (rx/2)-(rx/128), (ry/2)+(ry/128), (rx/2)+(rx/128) )		//if this results in problems it can also be set to fixed size (e.g. +/- 3 pixels)
		imagedisplay indisp = fft.imagegetimagedisplay(0)
		indisp.ImageDisplayAddROI(theROI)
		ROISetVolatile( theROI, 0 )
		ImageDisplaySetROISelected( inDisp, theROI, 1 )
		ROISetResizable( theROI, 0 )
		
		self.reset_UI()
	}

//apply a Hanning window to the input image or ROI, which reduces FT streaking
//inefficient memory useage since it creates a full copy of the stack
//for large stacks, this needs to be rewritten to apply the window to the active slice at each step
void hanning_filter(object self){
		image hann = slice3(front,l,t,0,0,rx,1,1,ry,1,2,sz,1)
		image hannrow := realimage("",4,rx,ry)*0
		image hanncol := realimage("",4,rx,ry)*0
		hannrow=0.5*(1-cos((2*pi()*(icol-(rx/2)))/(rx-1)))
		hanncol=0.5*(1-cos((2*pi()*(irow-(ry/2)))/(ry-1)))
		image hannimg = hannrow*hanncol
		hannimg -= min(hannimg)
		hannimg /= (max(hannimg)-min(hannimg))	//renormalization
		shiftcenter(hannimg)
		for(number i = 0; i<sz; i++){
			slice2(hann,0,0,i,0,rx,1,1,ry,1) *= hannimg
			}
		showimage(hann)		//push the new stack to the front display
		result("\nApplied Hanning window to original stack. Reload the original data to apply the final shifts!");
		log += "\nHanning window on real-space stack"
		self.load_stack()	//refresh active stack with the Hanning stack
	}

//aperture utility functions
//load an existing aperture image and apply it to the FT preview
void load_aperture(object self){	//load existing aperture into mask and apply to the ft
	try{
		string path
		if(!OpenDialog( path )) exit(0);
		mask = OpenImage(path)
		mask.hideimage()
		fft = oldfft*mask	//if the mask has the wrong dimensions this triggers the catch(), reverting to the previous mask
		posmask = mask*0
		result("\nMask image import from "+path+"")
		log += "\nLoaded existing mask "+getname(mask)+""	
		}
	catch{
		result("\nMask image import has wrong size for stack of size "+sx+"/"+sy+"")
		log += "\nAttempted to load existing mask "+getname(mask)+", wrong size"
		mask = oldmask
		}
	oldmask = mask
	}

//save current aperture to file
void save_aperture(object self){
	mask = tert(mask+posmask > 1, 1, mask)	//collapse the positive mask onto the image, and scale it to [0..1] again 
	showimage(mask)
	mask.ImageSetDescriptionText(log)	//also saving the log that led to this mask
	string path
	SaveAsDialog( "Save aperture to:", "FFT driftcorr mask ", path)
	SaveAsGatan( mask, path )
	mask.hideimage()
	log += "\nSaved current mask to "+path+""
	}

//circular centered top-hat aperture
//if the input image is not square, the mask is made elliptical and the radius field determines the y-size
void tophat(object self){
	hat_radius = DLGGetValue(self.lookupelement("hatradius"))
	number aspect = ry/rx		//rx/ry contain the ROI size or the full image in case of no ROI
	oldmask = mask
	mask = mask + tert(  sqrt( ((icol-rx/2)*aspect)**2 +  (irow-ry/2)**2 ) < hat_radius, 0, -1)
	mask = tert(mask < 0, 0, mask)		//some combinations (e.g. when a negative peak lies outside the radius) can result in values smaller 0, hence we cut off all negative values
	fft = oldfft*mask
	
	log += "\nApplied top-hat aperture with radius x/y "+(aspect*hat_radius)+"/"+(aspect*hat_radius)+""
	self.setElementIsEnabled("undobutton",1)
	}

//circular B'worth aperture with border smoothness of *order* (default = 4)
//originally from D.G. Mitchell (dmscripting.com)
void bworth_aperture(object self){
	hat_radius = DLGGetValue(self.lookupelement("hatradius"))
	image butterworthimg=realimage("",4,rx,ry)*0
	number bworthorder = 4
	number aspect = ry/rx	//elliptical aspect ratio in case of non-square image
	if(!GetNumber("B'worth order:",bworthorder,bworthorder)) exit(0);
	butterworthimg=1/(1+0.414*(  sqrt( ((icol-rx/2)*aspect)**2 +  (irow-ry/2)**2 ) /hat_radius)**(2*bworthorder))	//fixed halfpointconst=0.414 results in intensity=0.5 at *hatradius*
	butterworthimg = tert(butterworthimg < 5e-2,0,butterworthimg)	//sharp cutoff at the tail to prevent floating-point problems, 5e-2 is arbitrary
	
	oldmask = mask
	mask = tert(mask*butterworthimg > 1,1,mask*butterworthimg)	//stacking of mask and bworth including removing too large values
	fft = oldfft*mask
	
	log += "\nApplied B'worth aperture of order "+bworthorder+" with radius x/y "+(aspect*hat_radius)+"/"+(aspect*hat_radius)+""
	self.setElementIsEnabled("undobutton",1)
	}
	
//Gaussian blur of the aperture (mask) image by convolution with a 5x5 kernel
//choice of kernel size is abitrary
void gaussblur(object self){
	number stddev = 1.5
	if(!GetNumber( "Gaussian blur width:", stddev, stddev)) exit(0);
	
	image kernel=realimage("Gaussian kernel",4,5,5)
	kernel=1/(2*pi()*stddev**2)*exp(-1*(((icol-2)**2+(irow-2)**2)/(2*stddev**2)))	//irow-2 is for the case of a 5x5 kernel, generalized should be irow-sx/2-1 ? (dimensions start at 0)
	
	oldmask = mask
	mask = tert(mask+posmask > 1, 1, mask)	//collapse the positive spots onto the mask before blurring
	mask = convolution(mask,kernel)
	fft = oldfft*mask	
	
	log += "\nApplied Gaussian blur with standard deviation "+stddev+""
	self.setElementIsEnabled("undobutton",1)
	}

//dampen the mask at the position of a maximum within the FT ROI
void apply_negative_spot(object self){
	number lf,tf,rf,bf				//coordinates of the FT ROI
	if(!self.get_ROI_position(fft,lf,tf,rf,bf)){
		result("\nNo ROI in FFT preview image found, place a rectangle ROI around a spot of interest first!")	//in case the ROI was accidentally deleted
		exit(0)
		}
	number px,py
	max(fft[tf,lf,bf,rf],px,py)		//coarse maximum within the ROI
	try{ 
		self.imagerefineextrema(fft[tf,lf,bf,rf],px,py)	//sub-pixel refinement
		px += lf					//shift back the origin to the whole image
		py += tf
		px = max(px,lf)				//test whether the refined maximum is within the ROI, if not then shift it to the nearest ROI edge
		py = max(py,tf)
		px = min(px,rf)
		py = min(py,bf)
		if(verbose) result("\nPlacing mask around maximum at :"+px+"/"+py+"")
		oldmask = mask
		mask = self.pearson(mask,5,5,px,py)				//hard-coded width=5 seems to work well for centerband peaks
		fft = oldfft*mask
	
		log += "\nRemoved local maximum at "+px+","+py+""
		self.setElementIsEnabled("undobutton",1)
		}
	catch
		{
		if(verbose) result("\nFailed to place mask within the ROI, proceeding without any changes")	//this can happen if the FT ROI edge is placed exactly on the peak
		return
		}
	}

//positive circular mask collection
//these are hidden until the mask is saved
void apply_positive_spot(object self){
	number lf,tf,rf,bf				//coordinates of the FT ROI
	if(!self.get_ROI_position(fft,lf,tf,rf,bf)){
		result("\nNo ROI in FFT preview image found, place a rectangle ROI around a spot of interest first!")
		exit(0)
		}
	number px,py
	max(fft[tf,lf,bf,rf],px,py)		//coarse maximum within the ROI
	try{ 							//same logic as in the negative mask spot search
		self.imagerefineextrema(fft[tf,lf,bf,rf],px,py)	//sub-pixel refinement
		px += lf					//shift back the origin
		py += tf
		px = max(px,lf)				//test whether the refined maximum is within the ROI, if not then shift it to the nearest ROI edge
		py = max(py,tf)
		px = min(px,rf)
		py = min(py,bf)
		if(verbose) result("\nPlacing positive mask around maximum at :"+px+","+py+" within the ROI, this cannot be undone!")
		
		number spotradius = 5
		if(!GetNumber( "Positive mask radius:", spotradius, spotradius)) exit(0);
		posmask = tert(sqrt((icol-px)**2+(irow-py)**2) < spotradius,2,posmask)	//positive masks are circular (sharp) apertures around the maximum
																				//value of 2 guarantees that they collapse correctly onto the negative mask (in principle posmask = 1 and mask+posmask >= 1 would work as well)
		log += "\nAdded positive mask with radius "+spotradius+" at "+px+"/"+py+""
		//self.setElementIsEnabled("undobutton",1)	//positive spot masks cannot be undone (!)
		}
	catch
		{
		if(verbose) result("\nFailed to place mask within the ROI, proceeding without any changes")	//this can happen if the FT ROI edge is placed exactly on the peak
		return
		}
	}

//revert to previous state of the mask image
void undo_mask(object self){
	mask = oldmask
	fft = oldfft * mask
	
	log += "\nReverted last step"
	self.setElementIsEnabled("undobutton",0)	//can only undo one step
	}

//WIP functions specific to holography. Remove sideband center / remove Fresnel streak
//Note that the undo button also works for these compound masks
//Remove_sideband iterates over the full FFT and locates maxima within a certain search radius around the carrier frequency
void remove_sideband(object self){
	number px,py
	number width = 7		//sideband peaks are generally broader, so we increase the width of the negative peak
	
	hat_radius = DLGGetValue(self.lookupelement("hatradius"))	//get the hat radius as a sideband search window, optional
	oldmask = mask
	temp = oldfft
	temp = tert((irow < ry/2) && (iradius > 10),temp,0)			//top half of the modulus(FT) image minus the centerband pixel (hardcoded 10px iradius)
	max(temp,px,py)			//coarse maximum, in the case of an off-axis hologram with acceptable visibility this corresponds to the sideband carrier frequency
	if(verbose) result("\nSideband located (coarse) at :"+px+","+py+"");
	self.imagerefineextrema(temp,px,py)
	if(verbose) result("\nSideband located (fine) at :"+px+","+py+"");
	mask = self.pearson(mask,width,5,px,py)		//remove the sideband maximum with a larger width since it tends to be much more intense than the Bragg echoes
	mask = self.pearson(mask,width,5,rx-px,ry-py)
	temp = self.pearson(temp,width,5,px,py)		//applying the negative peak to the mask *and* the preview, the temp image is discarded later on	
	temp = self.pearson(temp,width,5,rx-px,ry-py)
	
	width = 5
	number counter = 4
	if(!GetNumber( "Number of sideband peaks to remove", counter, counter )) exit(0);
	log += "\nAttempted to remove "+counter+" local maxima within the sidebands"
	
	number search_radius = min( (sqrt((rx/2-px)**2+(ry/2-py)**2))*0.6 ,hat_radius )	//default guess = min(hat_radius field,0.6*carrier frequency)
	if(!getnumber("Sideband search radius:",search_radius,search_radius)) exit(0);
	if(verbose) result("\nSearch radius is : "+search_radius);
	temp = tert(sqrt((icol-px)**2+(irow-py)**2) < search_radius,temp,0)	//limit the sideband search range	
	
	while(counter > 0){
		max(temp,px,py)
		self.imagerefineextrema(temp,px,py)
		if(verbose) result("\nRemoving symmetric peak at :"+px+","+py+"")
		mask = self.pearson(mask,width,5,px,py)
		mask = self.pearson(mask,width,5,rx-px,ry-py)		//remove the symmetric peak on the -q_c side
		temp = self.pearson(temp,width,5,px,py)				//also remove the peak within the serach window... inefficient but otherwise mask and temp images might get mixed up
		temp = self.pearson(temp,width,5,rx-px,ry-py)
		counter -= 1
		}
		
	fft = oldfft*mask
	self.setElementIsEnabled("undobutton",1)
	}

//Remove_Fresnel is a line filter along the Fresnel streak, currently calling HW function
void remove_fresnel(object self){
		number px,py
		
		oldmask = mask
		temp = oldfft
		
		temp = tert((irow < ry/2) && (iradius > 10),temp,0)
		max(temp,px,py)	//sideband pixel position (absolute)
		self.imagerefineextrema(temp,px,py)
		if(verbose) result("\nLocated refined SB position at "+px+"/"+py+"");
		
		number Fresnel_length
		number Fresnel_width
		if(!getnumber("Relative length of Fresnel streak: ", 0.25, Fresnel_length)) exit(0);
		if(!getnumber("Pixel width of Fresnel streak: ", 7, Fresnel_width)) exit(0);
		px = px + (rx/2 - px)/2	//shifting the coordinates to the midpoint
		py = py + (ry/2 - py)/2
		
		number G1X,G1Y,G2X,G2Y	//start / endpoints of the streak, in absolute coordinates
		G1X = px + Fresnel_length * (rx/2 - px) * 2
		G1Y = py + Fresnel_length * (ry/2 - py) * 2
		G2X = px - Fresnel_length * (rx/2 - px) * 2
		G2Y = py - Fresnel_length * (ry/2 - py) * 2
		if(verbose) result("\nPlacing Fresnel filter from "+G1X+"/"+G1Y+" to "+G2X+"/"+G2Y+"");
		
		image Gauss_mask = HW_CreateGaussLineFilter( G1X, G1Y, G2X, G2Y, Fresnel_width, rx, ry )
		Gauss_mask *= Rotate(Gauss_mask,Pi())	//180deg rotation to mirror the +/- component
		
		mask *= Gauss_mask
		log += "\nPlaced Fresnel filter with "+Fresnel_length+" length, "+Fresnel_width+" width"
		
		fft = oldfft*mask
		self.setElementIsEnabled("undobutton",1)
	}

//UI button functions
//disable all buttons while thread is running
void block_UI(object self){

		self.setElementIsEnabled("getimagebutton",0)
		self.setElementIsEnabled("hanningbutton",0)
		self.setElementIsEnabled("apertureloadbutton",0)
		self.setElementIsEnabled("aperturesavebutton",0)
		self.setElementIsEnabled("tophatbutton",0)
		self.setElementIsEnabled("bworthbutton",0)
		self.setElementIsEnabled("posbutton",0)
		self.setElementIsEnabled("negbutton",0)
		self.setElementIsEnabled("blurbutton",0)
		self.setElementIsEnabled("undobutton",0)
		self.setElementIsEnabled("sidebandbutton",0)
		self.setElementIsEnabled("linefilterbutton",0)
		self.setElementIsEnabled("xcfbutton",0)
		self.setElementIsEnabled("pcfbutton",0)
		self.setElementIsEnabled("load_drift_button",0)
		self.setElementIsEnabled("apply_drift_button",0)
	}

//enable all buttons after first assignment or thread completion
void reset_UI(object self){

		self.setElementIsEnabled("getimagebutton",1)
		self.setElementIsEnabled("hanningbutton",1)
		self.setElementIsEnabled("apertureloadbutton",1)
		self.setElementIsEnabled("aperturesavebutton",1)
		self.setElementIsEnabled("tophatbutton",1)
		self.setElementIsEnabled("bworthbutton",1)
		self.setElementIsEnabled("posbutton",1)
		self.setElementIsEnabled("negbutton",1)
		self.setElementIsEnabled("blurbutton",1)
		self.setElementIsEnabled("undobutton",0)
		self.setElementIsEnabled("sidebandbutton",1)
		self.setElementIsEnabled("linefilterbutton",has_HW)		//Fresnel filter functionality currently attempts to call the line filter function from HW
		self.setElementIsEnabled("xcfbutton",1)
		self.setElementIsEnabled("pcfbutton",1)
		self.setElementIsEnabled("load_drift_button",1)
	}

//image alignment threads
//image alignment by cross-correlation of either subsequent pairs (pairs==1) or vs. the first slice (pairs==0)
void xcfthread(object self){
	image out := realimage("",4,2,sz)
	if(pairs) out.setname("Pairwise XCF drift sequence dy/dx");
	else out.setname("XCF Drift sequence dy/dx vs. first");
	outID = getimageID(out)
	
	number px,py
	image current, next
	image FT1,FT2, IFXcorr
	compleximage conjugate, Xcorr
	
	current := slice2(front,l,t,0,0,rx,1,1,ry,1)	//front slice, ROI or whole image
	FT1 = realFFT(current)
	
	for(number i = 0; i<sz-1; i++){
		OpenAndSetProgressWindow( "Aligning:", ""+(i+1) , "of "+(sz-1))
		if(pairs){ current := slice2(front,l,t,i,0,rx,1,1,ry,1); FT1 = realFFT(current); }
		next := slice2(front,l,t,i+1,0,rx,1,1,ry,1)
		//begin XCORR
		FT2 = realFFT(next)
		conjugate = complex(real(FT2), -1*imaginary(FT2))
		Xcorr = (FT1*conjugate)*mask				//apply the mask onto the cross-corr image
		IFXcorr = realIFFT(Xcorr)
		ShiftCenter(IFXcorr)
		IFXcorr = (IFXcorr-min(IFXcorr))/(max(IFXcorr)-min(IFXcorr))	//renormalizing
		//end XCORR
		max(IFXcorr[ry/2-maxy,rx/2-maxx,ry/2+maxy,rx/2+maxx],px,py)		//pixel max within the limiting centered search window on the IFFT; search window is square, might want to make it circular
		//result("\nInitial XCORR max (uncorrected): "+px+"/"+py)
		px += rx/2-maxx		//shift back the origin
		py += ry/2-maxy
		//result("\nInitial XCORR max : "+px+"/"+py)
		
		try{
			self.imagerefineextrema(IFXcorr,px,py)	//operate on the whole image to avoid edge problems, could be replaced by IFXcorr[ry/2-maxy,rx/2-maxx,ry/2+maxy,rx/2+maxx] again
			px -= rx/2		//shift the origin to the image center
			py -= ry/2
			}
		catch{
			px = out.getpixel(1,i)		//rather than have a "jumpy" drift series, we replace failed measurements with the previous values (linear drift approximation)
			py = out.getpixel(0,i)
			result("\nSlice "+(i+1)+" Xcorr failed.")
			break
			}
		
		if(verbose) result("\nRefined Xcorr vector "+(i+1)+" = "+px+","+py)
		out.setpixel(1,i+1,px)			//we write y/x not x/y, for historical reasons (interfacing with another script)
		out.setpixel(0,i+1,py)
		}
		
	showimage(out)
	
	log += "\nMeasured Xcorr series with pairwise = "+pairs+""
	out.ImageSetDescriptionText(log)
	self.reset_UI()
	self.setElementIsEnabled("apply_drift_button",1)
	return
	}

//identical to XCF but with PCF normalization
void pcfthread(object self){	
	image out := realimage("",4,2,sz)
	if(pairs) out.setname("Pairwise PCF drift sequence dy/dx");
	else out.setname("PCF Drift sequence dy/dx vs. first");
	outID = getimageID(out)
	
	number px,py
	image current, next
	image FT1,FT2, IFXcorr
	compleximage conjugate, Xcorr
	number scale
	image norm
	
	current := slice2(front,l,t,0,0,rx,1,1,ry,1)	//front slice, ROI or whole image
	FT1 = realFFT(current)
	
	for(number i = 0; i<sz-1; i++){
		OpenAndSetProgressWindow( "Aligning:", ""+(i+1) , "of "+(sz-1))
		if(pairs){ current := slice2(front,l,t,i,0,rx,1,1,ry,1); FT1 = realFFT(current); }
		next := slice2(front,l,t,i+1,0,rx,1,1,ry,1)
		FT2 = realFFT(next)
		
		scale = (1e-6)*(0.5*sum(modulus(FT1[ry/2-1,rx/2-1,ry/2,ry/2]))+0.5*sum(modulus(FT2[ry/2-1,rx/2-1,ry/2,rx/2])))	//that should be the correct average intenisty of the center pixel IF the center is unique i.e. rx+ry = even
		conjugate = complex(real(FT2), -1*imaginary(FT2))
		Xcorr = (FT1*conjugate)
		norm = modulus(Xcorr) + scale
		Xcorr = mask*(Xcorr/norm)
		IFXcorr = realIFFT(Xcorr)
		ShiftCenter(IFXcorr)
		max(IFXcorr[ry/2-maxy,rx/2-maxx,ry/2+maxy,rx/2+maxx],px,py)
		//result("\nInitial XCORR max (uncorrected): "+px+"/"+py)                                                          
		px += rx/2-maxx
		py += ry/2-maxy
		//result("\nInitial XCORR max : "+px+"/"+py)
		
		try{
			self.imagerefineextrema(IFXcorr,px,py)
			px -= rx/2
			py -= ry/2
			}
		catch{
			px = out.getpixel(1,i)
			py = out.getpixel(0,i)
			result("\nSlice "+(i+1)+" PCF failed.")
			break
			}
		
		if(verbose) result("\nRefined PCF vector "+(i+1)+" = "+px+","+py)
		out.setpixel(1,i+1,px)
		out.setpixel(0,i+1,py)
		}
	
	showimage(out)
	
	log += "\nMeasured PCF series with pairwise = "+pairs+""
	out.ImageSetDescriptionText(log)
	self.reset_UI()
	self.setElementIsEnabled("apply_drift_button",1)
	return
	}

//apply measured drift vectors via FT shift thread
void apply_thread(object self){
	image drift := GetimagefromID(outID)
	image out = imageclone(front)
	
	compleximage xgrid = compleximage("",8,sx,sy), ygrid = compleximage("",8,sx,sy)
	complexnumber j = complex(0,1)
	number offsetx = 0, offsety = 0, maxsx = 0, maxsy = 0
	image current
	compleximage acFFT, shift	//CFFT was being overloaded here
	
	xgrid = -j*2*pi()*(icol-sx/2)/sx
	ygrid = -j*2*pi()*(irow-sy/2)/sy
	
	for(number i = 0; i<sz; i++){
		OpenAndSetProgressWindow( "Applying:", ""+i+"" , "of "+(sz-1)+"")
		current = slice2(out,0,0,i,0,sx,1,1,sy,1)
		acFFT = realFFT(current)
		offsetx = (pairs*offsetx)+drift.GetPixel(1,i)*-1	//axes are flipped vs. Python, irrelevant for this script
		offsety = (pairs*offsety)+drift.GetPixel(0,i)*-1	//pairs is either 0 or 1, which results in us either accumulating values or not (working as intended)
		maxsx = max(abs(maxsx),abs(offsetx)) 				//calculating the maximum shift as symmetric, which removes more edge than necessary; could be done smarter
		maxsy = max(abs(maxsy),abs(offsety))				//set use_cut = 0 to do this manually
		//if(offsetx<0) maxsx *= -1							//if we want to make smarter choices about cutting the edges, we may want to store the sign						
		//if(offsety<0) maxsy *= -1
		result("\nShifting slice "+i+" by "+offsetx+"/"+offsety)
		//shift = exp( -j*2*pi()* ( ((xgrid*offsetx)/sx + (ygrid*offsety)/sy) ))	//no longer assuming sx=sy
		shift = exp(xgrid*offsetx+ygrid*offsety)	        //offset values are in absolute pixels, so independent of the grid size
		acFFT *= shift
		current = realIFFT(acFFT)
		slice2(out,0,0,i,0,sx,1,1,sy,1) = current
		}
	
	image out_cut 
	if(use_cut){
		out_cut = out[ceil(abs(maxsy)),ceil(abs(maxsx)),sy-ceil(abs(maxsy)),sx-ceil(abs(maxsx))]	//this can use some work
		out_cut.setname(getname(front)+"_aligned_cut")
		}
	else{
		out_cut = out
		out_cut.setname(getname(front)+"_aligned")
		}
	
	taggroup tgsrc, tgdest		//copy all metadata
	tgsrc = imagegettaggroup(front)
	tgdest = imagegettaggroup(out_cut)
	TagGroupCopyTagsFrom( tgdest, tgsrc )
	ImageCopyCalibrationFrom( out_cut, front )
	
	showimage(out_cut)
	log += "\nApplied drift series '"+getname(drift)+"' with pairwise = "+pairs+""
	out_cut.ImageSetDescriptionText(log)
	self.reset_UI()
	return
	}

//prepare the appropriate alignment threads and run
void crosscorr_thread(object self){
	maxx = floor(DLGGetValue(self.lookupelement("maxx")))
	maxy = floor(DLGGetValue(self.lookupelement("maxy")))
	pairs = DLGGetValue(self.lookupelement("pairs"))
	
	mask = tert(mask+posmask > 1, 1, mask)	//collapse the positive spots onto the mask
	
	self.block_UI()
	self.StartThread("xcfthread")
	}
	
void phasecorr_thread(object self){
	maxx = floor(DLGGetValue(self.lookupelement("maxx")))
	maxy = floor(DLGGetValue(self.lookupelement("maxy")))
	pairs = DLGGetValue(self.lookupelement("pairs"))
	
	mask = tert(mask+posmask > 1, 1, mask)	//collapse the positive spots onto the mask
	
	self.block_UI()
	self.StartThread("pcfthread")
	}

//load existing drift series, keep in mind the order (dy/dx along the y-dimension)
void load_drift(object self){
	image temp
	if(!getoneimagewithprompt("Drift sequence image:","Drift",temp)) exit(0)
	outID = getimageID(temp)
	self.setElementIsEnabled("apply_drift_button",1)
	}

//apply drift series thread call
void apply_drift(object self){
	pairs = DLGGetValue(self.lookupelement("pairs"))
	
	self.block_UI()
	self.StartThread("apply_thread")
	}


//build the interface
taggroup CreateMainDialog(object self){
		number buttonwidth = 44
		
		taggroup Dialog = DLGCreateDialog("Main Dialog")
		taggroup Dialogelements
		taggroup Dialoggroup = DLGCreateGroup(Dialogelements)
		taggroup button_items
		taggroup button_box = DLGCreateBox("Input",button_items)
		taggroup assignbutton = DLGCreatePushbutton("Assign","load_stack").DLGidentifier("getimagebutton").DLGwidth(buttonwidth)
		taggroup hanningbutton = DLGCreatePushbutton("Hanning","hanning_filter").DLGidentifier("hanningbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup buttongroup = DLGGroupItems(assignbutton, hanningbutton).dlgtablelayout(2,1,0)
		button_items.DLGAddElement(buttongroup)
		Dialog.DLGAddElement(button_box)
		
		taggroup aperture_items
		taggroup aperture_box = DLGCreateBox("Aperture",aperture_items)
		taggroup apertureloadbutton = DLGCreatePushbutton("Load","load_aperture").DLGidentifier("apertureloadbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup aperturesavebutton = DLGCreatePushbutton("Save","save_aperture").DLGidentifier("aperturesavebutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup aperturegroup = DLGGroupItems(apertureloadbutton, aperturesavebutton).dlgtablelayout(2,1,0)
		aperture_items.DLGAddElement(aperturegroup)	
		Dialog.DLGAddElement(aperture_box)
		
		taggroup mask_items
		taggroup mask_box = DLGCreateBox("Radial masks",mask_items)
		taggroup hatbutton = DLGCreatePushbutton("Top-hat","tophat").DLGidentifier("tophatbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup bworthbutton = DLGCreatePushbutton("B'worth","bworth_aperture").DLGidentifier("bworthbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup radial_mask_group = DLGGroupItems(hatbutton,bworthbutton).DLGtablelayout(2,1,1)
		mask_items.DLGAddElement(radial_mask_group)
		taggroup radiuslabel = DLGcreatelabel("Radius:")
		taggroup radiusfield = DLGCreateRealField(hat_radius,8,4).DLGIdentifier("hatradius")
		taggroup radius_group = DLGGroupItems(radiuslabel,radiusfield).DLGtablelayout(2,1,1)
		mask_items.DLGAddElement(radius_group)
		taggroup blurbutton = DLGCreatePushbutton("Blur","gaussblur").DLGidentifier("blurbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup undobutton = DLGCreatePushbutton("Undo","undo_mask").DLGidentifier("undobutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup positivebutton = DLGCreatePushbutton("Postive","apply_positive_spot").DLGidentifier("posbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup negativebutton = DLGCreatePushbutton("Negative","apply_negative_spot").DLGidentifier("negbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup misc_aperture_group = DLGGroupItems(positivebutton,negativebutton,blurbutton,undobutton).DLGtablelayout(2,2,1)
		mask_items.DLGAddElement(misc_aperture_group)	
		Dialog.DLGAddElement(mask_box)
		
		taggroup holo_items
		taggroup holo_box = DLGCreateBox("Holography",holo_items)
		taggroup sidebandbutton = DLGCreatePushbutton("Sideband","remove_sideband").DLGidentifier("sidebandbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup linefilterbutton = DLGCreatePushbutton("Fresnel","remove_fresnel").DLGidentifier("linefilterbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup hologroup = DLGGroupItems(sidebandbutton, linefilterbutton).dlgtablelayout(2,1,1)
		holo_items.DLGAddElement(hologroup)	
		Dialog.DLGAddElement(holo_box)
		
		taggroup thread_items
		taggroup thread_box = DLGCreateBox("Measure",thread_items)
		taggroup XCFbutton = DLGCreatePushbutton("XCF","crosscorr_thread").DLGidentifier("xcfbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup PCFbutton = DLGCreatePushbutton("PCF","phasecorr_thread").DLGidentifier("pcfbutton").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup pairswitch = DLGCreateCheckBox("Pairwise",1).DLGIdentifier("pairs").DLGenabled(1)
		taggroup maxxlabel = DLGcreatelabel("Max X")
		taggroup maxylabel = DLGcreatelabel("Max Y")
		taggroup xmaxfield = DLGCreateRealField(maxx,8,4).DLGIdentifier("maxx")
		taggroup ymaxfield = DLGCreateRealField(maxy,8,4).DLGIdentifier("maxy")
		taggroup maxgroup = DLGGroupItems(maxxlabel,maxylabel,xmaxfield,ymaxfield).dlgtablelayout(2,2,1)
		taggroup threadgroup = DLGGroupItems(xcfbutton, pcfbutton).dlgtablelayout(2,1,1)
		thread_items.DLGAddElement(threadgroup)
		thread_items.DLGAddElement(pairswitch)
		thread_items.DLGAddElement(maxgroup)
		Dialog.DLGAddElement(thread_box)
		
		taggroup drift_items
		taggroup drift_box = DLGCreateBox("Drift Sequence",drift_items)
		taggroup load_drift_button = DLGCreatePushbutton("Load","load_drift").DLGidentifier("load_drift_button").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup apply_drift_button = DLGCreatePushbutton("Apply","apply_drift").DLGidentifier("apply_drift_button").DLGenabled(0).DLGwidth(buttonwidth)
		taggroup driftgroup = DLGGroupItems(load_drift_button, apply_drift_button).dlgtablelayout(2,1,1)
		drift_items.DLGAddElement(driftgroup)	
		Dialog.DLGAddElement(drift_box)
		
		
		taggroup footer = DLGCreateLabel(version)
		Dialog.DLGAddElement(footer)
		
		DLGLayout(Dialog, DLGCreateTableLayout( 1, 7, 0 ))
		
		taggroup position			
		position = DLGBuildPositionFromApplication()
		position.TagGroupSetTagAsString( "Width", "Medium" )
		position.TagGroupSetTagAsString( "Side", "Top" )
		Dialog.DLGPosition(position) 

		return Dialog
	}


//constructor / deconstructor
DialogUI(object self){ self.super.init(self.CreateMainDialog()); }
~DialogUI(object self){ try deleteimage(fft); catch exit(0); }	//cleanup
}

void main()
{
	object dialog_UI = alloc(DialogUI)
	dialog_UI.init(dialog_UI.scriptobjectgetID())
	dialog_UI.display("Drift Corr UI")
}

main()
exit(0)