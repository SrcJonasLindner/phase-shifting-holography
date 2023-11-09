//An UI to aquire PS-EH tilt series using either the gun or beam tilt of the microscope
//
//Useage:
//
//"Start" reads the current field values and reponds by either acquiring a tilt series (using gun or beam tilt function via the "Gun" checkbox), or a static double-exposure series:
//
//if either the angle or amplitude stepsize fields are > 0, it records a tilt series with the appropriate ramp
//the "Pairs" checkbox determines whether the acquisition alternates between 0-tilt and current tilt step, or records a continuous series
//
//if both stepsize fields are set to 0, it records a static series of holograms at the given %startangle% / %startamp% values
//the %stopamp% field then determines the number of repetitions
//the "Pairs" checkbox again determines whether a series of 0-tilt + fixed-tilt pairs is recorded, or if only a series of holograms at a static tilt is recorded 
//
//in addition, the "Vacuum reference" checkbox determines whether the beam is shifted via UserImageShift to a reference position that has previously been defined with the "Ref." button
//this shift occurs in-between each step of the tilt series, the purpose is to enable parallel acquisition of a reference wave
//however, image shift comes with a small phase-shift of the hologram fringes, as well as some hysteresis; in practice, recording a sequential reference can be preferrable
//
//the "Calibrate" button allows manual input of calibration factors for the carthesian coordinate system of tilt-x and tilt-y
//it is advisable to measure these separately if you intend to use them; by default the amplitude fields scale the DAC values by 10e5, and single-digit integer values roughly correspond to micro-radians (at 300keV on a Titan G2)
//
//the "Wobble" button records a series of sawtooth-oscillations of the gun tilt, and the camera acquires the integral over multiple periods
//this requires entries in the %startamp%, %stopamp% and %stepamp% fields to determine the tilt amplitudes of the sawtooth oscillations that are sampled
//the %startangle% field determines the fixed direction, and the %stopangle% and %stepangle% fields do nothing (the %stepangle% variable is used to pass in the value in the "Samples" field, which determines the number of sampling points within one sawtooth)
//the end effect is that, for tilt amplitudes = n*2*pi() (where n is the samples field), each step in the sawtooth corresponds to a shift of one whole period, and the integrated pattern has maximum contrast
//in practice leaving the "Samples" field at 2.0 is sufficient and other values only require you to divide the amplitude by a different integer value to obtain the true 1*pi() phase shift
//see online documentation for further details
//
//"Measure stack" performs a measurement of the current front image stack, or of a square ROI in the image stack if one is present
//it attempts to perform a coarse drift measurement via cross-correlation, and also returns measured parameters (mean phase value and direction) of the sideband, as well as the standard deviation of intensities 
//in addition, it returns a 1D-image with two slices:
//the first slice contains the mean phase value of the sideband, which lets you estimate whether the biprism drift was sufficiently small / the tilt ramp properly chosen
//the second slice contains the standard deviation of intensities, which is used in the wobbler procedure to find the proper tilt amplitude for a phase shift of 2*pi() at the given biprism voltage, as described above
//(the standard deviation is also useful to evaluate the biprism stability / hologram visibility)
//
//if the parameter showprofile = 1 is set in the header of the script, the "Measure" button also returns an extracted profile series from the stack. the direction is chosen parallel to the sideband frequency
//for further use of these profiles, see the description of the "Refine" button Python integration in this script and the online repository
//
//"Drift" performs a phase-correlation function drift measurement on the current front image stack (T. Niermann, M. Lehmann, 2014, Micron 63, p. 28-34) with default parameters, and overwrites the drift values in the previous output (i.e. the coarse measurement needs to exist first)
//this functionality respects the "Pairs" checkbox, meaning that drift can be measured either vs. a fixed reference slice ("Pairs" off) or between subsequent slice pairs ("Pairs" on)
//"Correct" creates a copy of the current front image stack and applies the previously measured drift vectors
//
//"Refine" launches an external script file. What it does is up to you...
//
//Lastly, the driftx/drifty/phase fields show the average of the last-measured drift values and phase ramp 
//they are interactive since they used to be free parameters for linear drift correction, but currently only provide quick feedback on the quality of the measurement / evaluation
//
//***All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.
//***temscript is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENCE.txt file for any details.
//***All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.
//
//contact ross.ulrich at gmail.com for bug reports / feedback

//v0.10: added PCF drift measurement, image shift reference position 
//		 User Image Shift results in a small phase shift of the hologram fringes (!); thus calculating the reference wave from parallel acquisition still is not entirely correction-free
//		 PCF drift correction still has problems with sideband overlap, but is a decent first guess (much better than raw Xcorr)
//
//v0.14: fixed ROI adressing in measure_single
// v1.0: minor bugfixes, commenting and release
//
//known issues:	- the control interface for wobbler/tilt series/double exposure is more coplicated than it needs to be
//				- in order to measure PCF drift you need to highlight the image stack again (since each evaluation.init() grabs the current front image)
//				- PCF ROI cutting is currently removed, use the DriftCorr_UI script if you want to have full control over the drift correction process

//global parameters
number do_break = 0; number nil; number showprofile = 0;
number torad = pi()/180
string command_line= "C:\\Users\\supervisor\\Desktop\\BeamTiltApi\\\\GunTiltCmdLine.exe";				//path to the compiled SetGunTilt .exe
string beamTiltPath="C:\\Users\\supervisor\\Desktop\\BeamTiltApi\\ReadGunTiltCoord.exe";				//path to the compiled GetGunTilt .exe 
string GetImageShiftPath="C:\\Users\\supervisor\\Desktop\\BeamTiltApi\\ReadImageBeamShiftcmdLine.exe";	//path to the compiled GetUserImageShift .exe
string wobblerpath="C:\\Users\\supervisor\\Desktop\\Wobbler\\WobblerCmdLineSleep.exe";					//path to the compiled GunTiltWobbler .exe
string DMScriptPath="C:\\Users\\supervisor\\Desktop\\BeamTiltApi\\CallPythonScript_mod.s";				//path to an auxiliary script that will create the input for a Python script, which in turn will be executed in the Gatan Miniconda3 environment
																										//this is currently bound to the "Refine" button in the UI, what it does is up to you
																										//the example we provide (CallPythonScript_mod.s) passes an image ID and some other parameters to an input field in a Python script
																										//the example Python script (FitCosinusMain.py) uses the Gauss-Newton method to fit a cosine to a stack of profiles
																										//other uses one might imagine include performing the entire phase-shifting holography reconstruction, however this requires an aligned hologram series first
number SleepUntil = 50		//the wobbler functionality requires a delay time (ms) to avoid desynchronization between DM and the external process
string version = "v1.0 U. Ross, J. Lindner, Nov. 2023"

if(!getpersistentnumbernote("Tilt Series:LastID", nil))	//will be store for profile stack ID
	{
		setpersistentnumbernote("Tilt Series:LastID", 0)
	}

if(!getpersistentnumbernote("Tilt Series:OutID", nil))	//will be store for evaluation output ID
	{
		setpersistentnumbernote("Tilt Series:OutID", 0)
	}
	
setpersistentnumbernote("Tilt Series:ShiftX", 0)
setpersistentnumbernote("Tilt Series:ShiftY", 0)
//could also save the calibration values in global tags, but need to differentiate between gun and beam tilt calibration then

//thread communication interface
interface threadinterface
{
	void halt_acquisition(object self);
}

//the acquisition thread controls the gun/beam tilt wobbler as well as the hologram series acquisition
class acquisition_thread: Thread
{
number Dialog_ID

//parameters to be passed in from the UI
number double, guntilt, refswitch
number startangle,stopangle,stepangle
number startamp,stopamp,stepamp
number x0,y0,xcal,ycal
number shiftx, shifty

//communicate with the EM scripting API in order to read the current gun tilt
void getGunTilt(object self, number &x, number &y, number digits)
{
	number factor = 10**digits;		//significant digits of the unsigned long return value to be stored in a signed DM variable
	string command = beamTiltPath+ " " + factor + " 0" +"\n";	//fag 0 is x-coord
	x = LaunchExternalProcess( command )

	if (x > factor) {				//see documentation for GetGunTilt.exe for details
		x = (x-factor)*-1;			//first digit of the return value determins the sign
	}
	x /= factor;
	
	command = beamTiltPath+ " " + factor + " 1" +"\n";		//fag 1 is y-coord
	y = LaunchExternalProcess( command )
	if (y > factor) {
		y=(y-factor)*-1;
	}
	y /= factor;
}

//EM image shift setter, overwrites sx/sy with the initial values to undo the shift operation later
void apply_EMimage_shift(object self, number &sx, number &sy)
{
	number startX, startY
	EMGetImageShift( startX, startY )
	EMSetImageShift( sX, sY )
	result("\nImage shift set : "+sx+" / "+sy)	
	sx = startx; sy = starty;
}

//EM beam tilt setter, converts polar into carthesian coordinates and uses the internal DM function 
void apply_EMbeam_tilt(object self, number amp, number angle)
{
	number x = amp*cos(angle)
	number y = amp*sin(angle)
	EMSetBeamTilt(x0+(x/xcal), y0+(y/ycal))
}

//gun tilt setter, converts polar to carthesian and launches the external API process
void apply_gun_tilt(object self, number amp, number angle)
{
	number x = amp*cos(angle)
	x = (x0+(x/(xcal/0.00001))) //hard-coded scaling factor 10e-5 so the gun tilt values are approximately close to EMbeamtilt DAC values
	number y = amp*sin(angle)
	y = (y0+(y/(ycal/0.00001)))
	string command = command_line + " " + x + " " + y +"\n";
	number debug = LaunchExternalProcess(command);	//debug stores error codes, does not actually do anything currently
	result("\nGunTilt set : "+x+" / "+y)
}

//initialize the thread from the UI
//flag1 = single acquisition (0), use pairwise acquisition (current tilt + 0 tilt) (1) or launch the wobbler function instead (3)
//flag2 = gun tilt (1) or beam tilt (0)
//flag3 = shift to a reference position after each tilt step; usually not necessary since one can acquire an empty series afterwards; parallel acquisition minimizes biprism drift offset between reference and measurement but increases sample drift
//calx/y = calibration values; usually we operate without calibration, these have to be measured independently by e.g. the wobbler or by imaging the zero-disk shift in diffraction mode
//angle/amp fields are read from the UI, but these are treated slightly differently depending on context:
//
//	in the case of a series acquisition with %stepamp% > 0, this behaves as expected to acquire a tilt series
//
//	in the case of the gun tilt wobbler function call, the amplitude fields behave as expected (start, stop and stepsize of the wobbler amplitude series)
//	the %startangle% field determines the direction, and the %stepangle% field is used to forward the entry of the "wobbler samples" field define the sampling rate of the oscillator for each individual amplitude step 
//	(in the case of the gun tilt wobbler, the %stopangle% field does nothing)
//	see online documentation / publication for details on the wobbler calibration procedure
//
//	in case that both stepsize fields are 0, and a series acquisition is launched, a static series of tilt pairs (0 and current %startamp%/%startangle% field) is launched, and the %stopamp% field is used to determine the number of repetitions
//	this is an (inelegant) way to implement double-exposure holography without having to write a dedicated function
//
//	could probably be cleaned up in future versions to make the interface more intuitive
object init(object self, number flag1, number flag2, number flag3, number calx, number caly, number anglestart, number anglestop, number anglestep, number ampstart, number ampstop, number ampstep)
{
	double = flag1; guntilt = flag2; refswitch = flag3;
	startangle = anglestart*torad; stopangle = anglestop*torad; stepangle = anglestep*torad;
	startamp = ampstart; stopamp = ampstop; stepamp = ampstep;
	xcal = calx; ycal = caly;
	if(guntilt){
		result("\nUsing gun tilt.")
		self.getGunTilt(x0,y0,7)	//read the current gun tilt value as reference 0 tilt
		}
	else{
		result("\nUsing bound EMBeamTilt method.")
		EMGetBeamTilt(x0, y0)		//read the current beam tilt value as reference 0 tilt
		}
	result("\nZero tilt setting read : "+x0+" / "+y0)
	if(refswitch){					//read the previously stored reference image shift position
		getpersistentnumbernote("Tilt Series:ShiftX", shiftX)
		getpersistentnumbernote("Tilt Series:ShiftY", shiftY)
		}
	return self
}

//link the thread ID to the main UI thread
void link_thread(object self, number threadID)
{
	Dialog_ID = threadID
}

//prepare the camera for acquisition
//reads the current parameters for "acquire" in the DM user interface and returns a parameter set object
object prepare_cam(object self, object &cam, number &sx, number &sy)
{
	number binx, biny
	cam = CM_GetCurrentCamera()
	object acq_params=CM_GetCameraAcquisitionParameterset(cam, "Imaging", "Acquire", "Record", 1)
	CM_Validate_AcquisitionParameters(cam, acq_params)
	//CM_SaveCameraAcquisitionParameterSet(cam, acq_params) 	//updates the parameter set in the global tags, buggy
	CM_CCD_GetSize(cam,sx,sy)
	CM_GetBinning( acq_params, binx, biny )
	sx = sx/binx; sy = sy/biny;		//recalculate the size to the binned size so that the output size matches; this is passed back by reference (&sx/&sy)
	return acq_params
}

//gun wobbler external process call & sync acquisition
void gunwobbler(object self)
{
	object cam
	number sx,sy,sz
	string name
	object called_dialog = GetScriptObjectFromID(Dialog_ID) //the Dialog_UI object in case the UI is used to interrupt the acquisition
	
	//prepare the passed in variables for the gun wobbler case
	stepangle = round(stepangle/torad)								//we use %stepangle% to set the number of samples (periods) of the sawtooth oscillation, thus need to bring it back to integer value
	startamp *= 0.00001; stopamp *= 0.00001; stepamp *= 0.00001; 	//scale amp values by 10e-5 to approximately match EMBeamTilt DAC values
	if(startamp == 0) startamp = 0.000001 							//workaround until API catches 0
	number iters = 1+trunc((stopamp-startamp)/stepamp)				//number of sampling points in the amplitude series
	
	//determine the actual sawtooth time via API
	string command = wobblerpath+ " " + stepangle + " " +startangle +" "+ (startamp+stepamp) + " " + 10+ " " + SleepUntil; //using a default number of cycles of 10 to measure the mean sawtooth time
	number meanSawToothTime = LaunchExternalProcess(command);
	number meansampletime = meansawtoothtime/(stepangle*2)
	result("\nMean saw tooth period [ms]: " +  meanSawToothTime + "\n")
	result("\nMean time per step [ms]: " +  meansampletime + "\n")
	self.apply_gun_tilt(0,0)	//reset gun tilt in case it hasnt happened in the API call
	
	//prepare the acquisition parameters and determine the output dimensions
	object acq_params = self.prepare_cam(cam,sx,sy)
	
	//determine optimum exposure in [s] from sawtooth time in [ms]
	number exposure = CM_GetExposure(acq_params)*1000
	exposure = (exposure - mod(exposure,meanSawToothtime))/1000	//adjust the exposure to more closely match the cycle time
	result("\nOptimum exposure = "+exposure)
	
	//increment num_cycles until > exposure so that it is guaranteed the wobbler doesnt terminate before acquistion + overhead is finished 
	number num_cycles = 0
	try{
		while((num_cycles*meanSawToothTime/1000)<(exposure*2)){ 	//exposure*2 is arbitrary, could be made smaller as long as it is > exposure + sleep time; is used to guarantee that CM_acquire() does not finish after LaunchExternalProcessAsync()
		num_cycles++
		}
	}
	catch{ result("\nError when calculating required number of cycles"); num_cycles = 5; return; }	//default to 5 cycles, should never happen
	
	//start acquisition
	image out := realimage("",4,sx,sy,iters+2) //check sz, might be 1 off
	out.setname("Wobbled acq: amp "+startamp + "to"+ stopamp+" angle "+(startangle/torad)+" sleep delay "+SleepUntil+" cycle time "+(meanSawToothTime*2)+" sampling rate "+stepangle)
	out.ImageSetDescriptionText("\Tilt wobbler series acquistion:\nStartamp: "+startamp+"\nStopamp: "+stopamp+"\nStepamp: "+stepamp+"\nAngle: "+startangle+"\nSampling points: "+stepangle+"\nRate: "+meanSampleTime+"\Cycles: "+num_cycles)
	out.showimage()
	image firstimg := CM_AcquireImage(cam,acq_params) 	//grab a 0-tilt-before reference, write it into the last slice to be consistent with the other acquistion functions
	slice2(out,0,0,iters+1,0,sx,1,1,sy,1) = firstimg
	out.ImageCopyCalibrationFrom(firstimg)				//transfer metadata
	TagGroupCopyTagsFrom(imagegettaggroup(out),imagegettaggroup(firstimg))
	
	
	number amp, processID
	
	for (number j=0; j<=iters; j++ ){
	
		if(do_break){ self.apply_gun_tilt(0,0); called_dialog.halt_acquisition(); } 	//the break condition to respond to the stop button is checked before each iteration 
																						//gun tilt *should* already be reset, but repeat it to be on the safe side
		amp = startamp+(j)*stepamp
		command = wobblerpath+ " " + stepangle + " " +startangle +" "+ amp + " " + num_cycles+ " " + SleepUntil;
		processID = LaunchExternalProcessAsync(command)
		sleep(0.1)																		//creates some artificial overhead for the process to start properly, might not be necessary
		//result("\nAsync process passed control " + j + "/" +iters+"\n" )
		slice2(out,0,0,j,0,sx,1,1,sy,1) = CM_AcquireImage(cam,acq_params)
		OpenAndSetProgressWindow( "Acquisition:", j+" of ", ""+(iters) )
		WaitForExternalProcess(processID)												//due to the num_cycles calculation, the process lives at least 2*exposure [ms]
		self.apply_gun_tilt(0,0) 														//reset gun tilt after each loop iteration
	}
	setpersistentnumbernote("Tilt Series:OutID", GetImageID(out))						//store the outputID in global tags
}

//gun or beam tilt discriminator
void tilt_gun_or_beam(object self, number amp, number angle)
{
	if(guntilt==1){
		self.apply_gun_tilt(amp, angle)
	}
	else{
		self.apply_EMbeam_tilt(amp, angle)
	}
}

//the StartThread response to launch the whole thing
void RunThread(object self)
{
	object called_dialog=GetScriptObjectFromID(Dialog_ID)
	if(called_dialog.ScriptObjectIsValid())
	{
		if( double == 3){ self.gunwobbler(); called_dialog.halt_acquisition(); } //starting and finishing the wobbler acquisition (by setting do_break = 1) here to be sure that RunThread exits to be called from the DialogUI thread
																				 //instead, it would be sensible to launch the acquisition thread with StartThread(method) instead, but this would require a minor rewrite
		object cam
		number sx,sy,sz
		number switch	//case discriminator for acquisition depending on input configuration
		string name
		number i,j
		number refshiftx = shiftx; number refshifty = shifty;
		
		object acq_params = self.prepare_cam(cam,sx,sy)
		//sx = 512; sy = 512; // !!!!!during offline testing the image size and binning does not match the dummy output, check the fixed size and adjust accordingly!!!!! 
		
		//verify the image sizes sz to avoid empty slices / out-of-bounds errors; very prone to bugs
		if(stepangle != 0 && stepamp == 0){ sz = (double+1+refswitch)*floor((stopangle-startangle)/(stepangle))+2; switch = 0; i=startangle; name = "Tilt series from "+(startangle/torad)+" to "+(stopangle/torad); }	//loop over angles
		if(stepangle == 0 && stepamp != 0){ sz = (double+1+refswitch)*floor((stopamp-startamp)/(stepamp))+2; switch = 1; i=startamp; name = "Tilt series from "+startamp+" to "+stopamp; }								//loop over amps
		if(stepangle == 0 && stepamp == 0){ stopamp-=1; sz = (double+1+refswitch)*(stopamp)+2; switch = 2; i=0; name = "Tilted to "+startamp+" (amp) / "+startangle+" (angle) and back, "+stopamp+" times"; } 			//static amp/angle pair with 0 for double-exposure
		
		image out := realimage("",4,sx,sy,sz)
		out.setname(name)
		out.showimage()
		j = 0
		
		while(!do_break)	//check manual interrupt condition at the start of each loop
		{
			try
			{
				OpenAndSetProgressWindow( "Acquiring:", j+" of ", ""+(sz-1) )
				if(switch == 0)	//iterating over the angle only
				{
					self.tilt_gun_or_beam(startamp,i)
					i = i+stepangle
					if(i>stopangle) do_break = 1
				}
				if(switch == 1)	//iterating over the amp only
				{
					self.tilt_gun_or_beam(i,startangle)
					i = i+stepamp
					if(i>stopamp) do_break = 1
				}
				if(switch == 2)	//only one static amp/angle
				{
					self.tilt_gun_or_beam(startamp,startangle)
					i++
					if(i>stopamp) do_break = 1 //%stopamp% determines the number of static iterations, as described above
				}
				slice2(out,0,0,j,0,sx,1,1,sy,1) = CM_AcquireImage(cam,acq_params)
				OpenAndSetProgressWindow( "Acquisition:", j+" of ", ""+(sz-1) )
				j++
				
				if(double){		//0-tilt image for double acquisition
					self.tilt_gun_or_beam(0,0)
					slice2(out,0,0,j,0,sx,1,1,sy,1) = CM_AcquireImage(cam,acq_params)
					j++
					}
					
				if(refswitch){	//shifting to reference position and back
					self.apply_EMimage_shift(refshiftx,refshifty)
					slice2(out,0,0,j,0,sx,1,1,sy,1) = CM_AcquireImage(cam,acq_params)
					self.apply_EMimage_shift(refshiftx,refshifty)
					refshiftx = shiftx; refshifty = shifty;
					j++
					}							
			}
			catch	//remove the try/catch to view the error message if acquistion fails
			{
				result("\nException in acquistion thread")
				if(j > sz) result("\nAttempted to acquire more images than available slices") 
				self.tilt_gun_or_beam(0,0)
				called_dialog.halt_acquisition()
				break
			}
		}
		
		self.tilt_gun_or_beam(0,0)	//tilt back to 0
		if(j>0){					//if the acquisition did not run into an early exception record one final 0-tilt image and copy the metadata
			image lastimg := CM_AcquireImage(cam,acq_params) 
			slice2(out,0,0,sz-1,0,sx,1,1,sy,1) = lastimg	//verify sz
			OpenAndSetProgressWindow( "Acquisition:", j+" of ", ""+(sz-1) )
			out.ImageCopyCalibrationFrom(lastimg)
			TagGroupCopyTagsFrom(imagegettaggroup(out),imagegettaggroup(lastimg))
			//result("\nCopied tags")
		}
		out.ImageSetDescriptionText("\Tilt series acquistion:\nStartamp: "+startamp+"\nStopamp: "+stopamp+"\nStepamp: "+stepamp+"\nStartangle: "+startangle+"\nStopangle: "+stopangle+"\Stepangle "+stepangle+"\nPairwise: "+double)	//store the parameters in the image info text
		out.showimage()
		//setpersistentnumbernote("Tilt Series:OutID", GetImageID(out))	//outID gets used later on for the evaluation output, so i'd rather deal with getfrontimage() in the evaluation thread init
		result("\nTilt aquisition done.")
		called_dialog.halt_acquisition()
	}
}

//constructor/destructor
acquisition_thread(object self)
{
}

~acquisition_thread(object self)
{
}

}

//the evaluation thread controls the various post-processing functions
class evaluation_thread: Thread
{
number Dialog_ID

image front
number switch 					//0 = stack, 1 = pairwise stack, 2 = drift measurement, 3 = correct drift by last measured values in outputID
number driftx,drifty,driftphase	//return values to write the UI fields, or input values to perform linear drift correction (outdated)
number t,l,b,r					//ROI vertices 0 and 2, or whole image size
number SBAngle, SBLength		//sideband parameters
number rx, ry					//centerpoint of the extracted profile


//function to return the spanning coordinates of the first ROI, or of the whole image
number get_ROI_position(object self, image in, number &px, number &py, number &sx, number &sy)
{
	imagedisplay imgdisp=imagegetimagedisplay(in,0)
	
	try	
	{
		ROI theROI = ImageDisplayGetROI( imgDisp, 0 )
		ROIisvalid(theROI)
		ROIGetVertex( theROI, 0, px, py )
		ROIGetVertex( theROI, 2, sx, sy )
	}
	catch
	{
		getsize(in,sx,sy)
		px = 0; py = 0;
		result("\nROI invalid, returning whole image")
		return(0)
	}
	result("\nROI found at t/l/b/r: "+py+"/"+py+"/"+sy+"/"+sx)
	return(1)	
}

//initialize the thread from the UI
//the current image stack is updated by getfrontimage() in the DialogUI thread before passing it forward here
//flag1 = measure metrics (drift/phase etc) vs. first (0), measure metrics pairwise (1), measure PCF drift(2), correct PCF drift (3)
//xdrift/ydrift/phasedrift and read from and written to the UI fields, but only give a rough estimate for mean linear drift in x/y direction and phase ramp 
object init(object self, image in, number flag1, number xdrift, number ydrift, number phasedrift)
{
	front := in
	switch = flag1
	driftx = xdrift; drifty = ydrift; driftphase = phasedrift;
	self.get_ROI_position(front,l,t,r,b) 		//get the position of a rectangular ROI in the source stack
												//using a ROI accelerates the measurements, and for the determination of the phase ramp a ROI should always be chosen in the vacuum region
												//ROIs do not have to be square, but some functions involving the FFT may not work as expected for non-square ROIs
	SBAngle = 0; SBLength = 0; rx = 0; ry = 0;	//reset the sideband and profile parameters
	return self
}

//link the thread ID to the main UI thread
void link_thread(object self, number threadID)
{
	Dialog_ID = threadID
}

//calculate standard deviation of image intensities; could be written in two lines but more readable here
number stabw(object self, image in)
{
	number mean
	number sx,sy; getsize(in,sx,sy)	//make sure that the input is 2D
	mean = mean(in)
	image temp = (in-mean)**2
	return sqrt(sum(temp)/(sx*sy))
}

//cross-wise parabolic fit, originally by D.G. Mitchell (dmscripting.com) and adjusted for our purposes
number ImageRefineExtrema(object self, image img, number &px, number &py)
{
		number tx1, tx2, ty1, ty2
		number numa, numb, numc
		number valX, valY

		// do a parabolic fit in direction x
		numa = img.GetPixel(px-1, py)
		numb = img.GetPixel(px  , py)
		numc = img.GetPixel(px+1, py)
		tx2=(numa-numc)/2.0
    	tx1=(numa+numc)/2.0-numb
		valX=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb

		// do a parabolic fit in direction y
		numa = img.GetPixel(px, py-1)
		numb = img.GetPixel(px, py  )
		numc = img.GetPixel(px, py+1)
		ty2=(numa-numc)/2.0
    	ty1=(numa+numc)/2.0-numb
		valY=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb

		// update pixel position
		px += tx2/(2.0*tx1)
		py += ty2/(2.0*ty1)
		
		// check to see whether to return mimimum or maximum
		If(numa>=numb && numc>=numb) return min(min(valX, valY), numb)
		If(numa<=numb && numc<=numb) return max(max(valX, valY), numb)

		Return numb
}

//cross-correlation for coarse drift measurement
image CrossCorrelateWithSubPixelPrecision(object self, image src, image ref, number &pX, number &pY, number &qval)
{
		Number sizeX,sizeY
		src.GetSize(sizeX,sizeY)

		RealImage xcorr = src.CrossCorrelate(ref)
		
		// find the maximum and its coordinates in correlation image		
		Number spotX, spotY
		//qval = max(xcorr[sizeY/4,sizeX/4,3*sizeY/4,3*sizeX/4], spotX, spotY)	//limiting the maximum search to the passed in image center, have to adjust the pixel shift to center if used
		qval = max(xcorr, spotX, spotY)
	
		// refine with sub-pixel accuracy
		try(qval = self.ImageRefineExtrema(xcorr, spotX, spotY))
			catch {
				spotY = 0	//default to 0 (image center) if the maximum shifts out of the image frame 
				spotY = 0
				qval = 0
				result("\nCross-correlation center out of bounds")
				break
				}
		pX = spotX-(sizeX/2) //-sizeX/4 //shift pixel origin to image center
		pY = spotY-(sizeY/2) //-sizeY/4

		return xcorr	//returns the cross-correlation image, which is typically discarded
}

//find sideband pixel position (not frequency) relative to image center (FT origin), takes an FFT image as input (typecasting in = compleximage)
number Find_sideband(object self, compleximage in, number &px, number &py)
{
	number sx,sy
	in.GetSize(sx,sy)
	realimage tertimg = tert(iradius>2 && irow<(sy/2),modulus(in),0) 				//cutting off the bottom part of the FFT and the center band maximum (+2 pixels radius) to make the sideband position unique
	max(tertimg,px,py)																//for very low visibility, the sideband peak may not be the brightest peak in the modulus image; but in that case the entire hologram has to be discarded
	number pval = phase(getpixel(in,px,py))
	try{ self.imagerefineextrema(tertimg,px,py); px = (px-sx/2); py = (py-sy/2); }	//subpixel sideband position relative to FFT center
	catch{ result("\nRefining sideband position failed."); px = (px-sx/2); py = (py-sy/2); }		//if px/py = 0 the subsequent code may throw NAN/INF, which is less than ideal; we return the coarse pixel position relative to the image center, but this may still be erroneous
	return pval				//returns the coarse phase value; since the complex sideband pixel is a pole, interpolating the phase value wouldn't make sense; thus the example of the cosine refinement in real-space via Python
}

//measure a single slice in the series
//
//returns an extracted line profile for further processing (e.g. with the Python hook)
//this profile is centered on the ROI or the whole image, and the length is set to the minimum of ROI diameter or 9*sideband frequency
//the integration width his hard-coded to 32 pixels (regardless of ROI size)
//
//the metrics returned are:
//offsetx/y = Xcorr shift
//offsetval = Xcorr value	(later on discarded)
//pval = coarse phase of the sideband
//px/py = sideband vector in frequency space
//stab = standard deviation of intensities
//
// !!! possible bug in the calculation of the sidband vector/length !!!
// !!! ROI profiles may reach out of the image dimensions if the ROI is too close to the image edge !!!
// not critical since the extracted profile / tx/ty values for PS-EH reconstruction are advanced features that should be evaluated separately as well
image measure_single(object self, image in1, image in2, number &offsetx, number &offsety, number &offsetval, number &pval, number &px, number &py, number &stab) 
{
	number sx,sy
	in1.getsize(sx,sy)
	compleximage FFTimage = realFFT(in2[t,l,b,r])	//in2 is measurement slice, in1 is reference slice
	self.CrossCorrelateWithSubpixelPrecision(in2[t,l,b,r],in1[t,l,b,r],offsetx,offsety,offsetval) //throwing away the crosscorr image; now using ROI images
	pval = self.Find_sideband(FFTimage,px,py); px = px/(r-l); py = py/(b-t); //get sideband pixels, convert to frequency
	stab = self.stabw(in2[t,l,b,r])
	
	//finding the sideband angle  and length if it hasnt been initialized yet
	number anglemax, anglemin, angle, length
	if(!SBAngle){
		anglemax = max(py,px)
		anglemin = min(py,px)
		angle = atan2(anglemin/anglemax,anglemax/anglemax) //verify that the order of operations is correct
		result("\nSideband reference angle : "+(angle/torad)+" [deg]")
		SBangle = angle
	}
	else{
		angle = SBAngle
	}
	if(!SBLength){
		SBlength = (1/(sqrt(px**2+py**2))) 				//hardcoding profile length of 9*wavelength, px/py are sideband vector components; not 100% sure this calculation is correct, possibly of by a factor of 2
		length = 9*SBlength
		result("\nPx/Py = "+px+" , "+py+", pi() = "+SBlength+" pixels")
	}
	else{
		length = 9*SBLength
	}
	
	if(!rx || !ry){		//if profile has not been initialized; resets on self.init()
		if(r-l == sx && b-t == sy){ rx = sx/2; ry = sy/2; }	//if no ROI present, place the profile vertices around the image center
		else{
			rx = l+(r-l)/2; ry = t+(b-t)/2	//center pixel of the ROI in the front slice
		}
		//result("\nProfile length : "+length+" , Profile center: "+rx+","+ry)
	}
	
	number exs,exe,eys,eye
	if(length>(r-l) || length>(b-t)) length = min((r-l),(b-t))  //guarantee that length is not larger than the smallest circle in the ROI
	exs = rx+cos(angle)*length/2								//calculate the line profile vertices
	eys = ry+sin(angle)*length/2
	exe = rx-cos(angle)*length/2
	eye = ry-sin(angle)*length/2
	image profileimg = LiveProfile_ExtractLineProfile( in2, exs, eys, exe, eye, 32 )	//32 pixel integration width hard-coded
																						//a profile of 32px width will still run out of the ROI (although it is hard to predict by how much depending on the angle and ROI radius)
																						//make sure not to place the ROI immediately at the image edge if this causes problems, or shorten the profile
	return profileimg
}

//write measure_single results to an existing output image
void write_results(object self, image &out, number i, number offsetx, number offsety, number pval, number px, number py, number stddev)
{
	try{
	out.setpixel(i,0,offsetx)
	out.setpixel(i,1,offsety)
	out.setpixel(i,2,pval)
	out.setpixel(i,3,px)
	out.setpixel(i,4,py)
	out.setpixel(i,5,stddev)
	}
	catch{
	result("\nError writing results, probably output shape not matching passed-in parameters.")	//this shouldn't ever happen, but if it does the index i is probably wrong
	return
	}
}

//Fourier-shift image by drift vector
image FFT_shift(object self, image in, number px, number py)
{
	number sx,sy
	in.GetSize(sx,sy)
	image xgrid = realimage("",4,sx,sy), ygrid = realimage("",4,sx,sy)
	xgrid = icol-sx/2
	ygrid = irow-sy/2
	complexnumber j = complex(0,1)
	
	compleximage cFFT = realFFT(in)
	compleximage shift = exp( -j*2*pi()* ( ((xgrid*px)/sx + (ygrid*py)/sy) ))	// no longer assuming sx=sy
	cFFT *= shift
	realimage current = realIFFT(cFFT)
	return current
}

//construct FFT aperture, see D.G. Mitchell's website
image Bworth_aperture(object self, image in, number radius, number bworthorder)
{
	number sx,sy
	in.GetSize(sx,sy)
	image butterworthimg=realimage("",4,sx, sy)
	butterworthimg=0
	number halfpointconst=0.414

	butterworthimg=1/(1+halfpointconst*(iradius/radius)**(2*bworthorder))
	return butterworthimg 
}

//find drift by normalized Xcorr (PCF)
//for PCF reference see: T. Niermann, M. Lehmann, 2014, Micron 63, p. 28-34
image Xcorr_norm(object self, image src, image dest, number &px, number &py, image aperture)
{
	number sx,sy
	src.GetSize(sx,sy)
	image FTs, FTd
	FTs = realFFT(src)
	FTd = realFFT(dest)
	compleximage conjugate = complex(real(FTd),(-1)*imaginary(FTd))
	compleximage Xcorr = FTs*conjugate
	number scale = (1e-6)*(0.5*sum(modulus(FTs[sx/2-1,sy/2-1,sx/2,sy/2]))+0.5*sum(modulus(FTd[sx/2-1,sy/2-1,sx/2,sy/2])))	//hard-coded epsilon = 1e-6; mean of centerband DC amplitude (source + dest)
																															//may require even dimensions, adjust if you want the more general case
	realimage norm = modulus(Xcorr) + scale

	Xcorr = aperture*(Xcorr/norm) //renormalization
	image out = realIFFT(Xcorr)
	
	ShiftCenter(out) //quadrant swap
	max(out,px,py)
	//val = IG12.IUImageFindMax( sy/4, sx/4, 3*sy/4, 3*sx/4, px, py, 1) //alternative subpixel max finder
	try{ self.imagerefineextrema(out,px,py); px = (px-sx/2); py = (py-sy/2); return out; }	//subpixel drift position relative to Xcorr center
	catch{ result("\nRefining drift vector failed."); px = 0; py = 0; out = 0; return out; }	//defaulting back to zero shift if the refinement fails
}

//the StartThread response to launch the various evaluation procedures
//similar to before, this could be split into different methods instead of the if(flag) then() useage, but would require a rewrite of the code
void RunThread(object self)
{
	object called_dialog=GetScriptObjectFromID(Dialog_ID)
	if(called_dialog.ScriptObjectIsValid())		//check that the UI is still open; shouldn't ever really fail since the hand-over is fast
	{	
		image first, current, profile
		number offsetx,offsety,offsetval,pval,px,py,stddev
		number sx,sy,sz

		if(switch == 0 || switch == 1){	//case 0 and 1 are measure stack with pairwise off/on
		try{			
			//prepare the output images
			result("\nEvaluation output shape: slice number :: driftx, drifty, phase value, sideband frequency vector x/y, std. dev., refined phase value")	//refined phase value is a work in progress and not filled by the quick evaluation	
			get3dsize(front,sx,sy,sz)
			image out := realimage("Evaluation results (raw)",4,sz,7)

			first = slice2(front,0,0,0,0,sx,1,1,sy,1)	//take any slice to prepare the profile stack, in this case the first slice (tilt0)
			profile = self.measure_single(first,first,offsetx,offsety,offsetval,pval,px,py,stddev)	//prepare one profile image to get the dimensions
			number spx,spy; getsize(profile,spx,spy)
			image profilestack := realimage("",4,spx,spy,sz)
			
			if(switch == 0){	//if not pairwise, return measurement vs. first
				
				for(number i = 0; i<sz; i++){
					current = slice2(front,0,0,i,0,sx,1,1,sy,1)		//so far we still operate on the whole image, only the measure_single function cuts down the image to the ROI; might gain some speed by only passing the ROI, but this already works fast enough
					OpenAndSetProgressWindow( "Working on:", (i+1)+" of ", ""+sz )
					profile = self.measure_single(first,current,offsetx,offsety,offsetval,pval,px,py,stddev)
					slice2(profilestack,0,0,i,0,spx,1,1,spy,1) = profile
					self.write_results(out,i,offsetx,offsety,pval,px,py,stddev)
				}
				
				DLGvalue(called_dialog.lookupelement("driftx"),sum(out[0,0,1,sz])/sz)	//updating the UI fields with the mean values
				DLGvalue(called_dialog.lookupelement("drifty"),sum(out[1,0,2,sz])/sz)
				DLGvalue(called_dialog.lookupelement("driftphase"),sum(out[2,0,3,sz])/sz)
			}
			
			if(switch == 1){	//if pairwise, evaluate drift between measurement slices (i) and reference slices (i+1)
				
				for(number i = 0; i<(sz-1); i+=2){					//iterating over all the measurement slices, measuring relative to next (tilt back 0)
					first = slice2(front,0,0,i+1,0,sx,1,1,sy,1) 	//current tiltback slice
					current = slice2(front,0,0,i,0,sx,1,1,sy,1)		//current measurement slice
					OpenAndSetProgressWindow( "Working on:", ((i/2)+1)+" of ", ""+sz )
					profile = self.measure_single(first,current,offsetx,offsety,offsetval,pval,px,py,stddev)
					slice2(profilestack,0,0,i,0,spx,1,1,spy,1) = profile
					self.write_results(out,i,offsetx,offsety,pval,px,py,stddev)
				}
				
				first = slice2(front,0,0,sz-1,0,sx,1,1,sy,1)		//last tiltback slice
				
				for(number i = 1; i<(sz-1); i=i+2){					//iterating over all the reference slices, measuring relative to the last reference slice
					current = slice2(front,0,0,i,0,sx,1,1,sy,1)		//current tiltback slice
					OpenAndSetProgressWindow( "Working on:", (((sz-1)/2)+((i/2)+1))+" of ", ""+sz )
					profile = self.measure_single(first,current,offsetx,offsety,offsetval,pval,px,py,stddev)
					slice2(profilestack,0,0,i,0,spx,1,1,spy,1) = profile
					self.write_results(out,i,offsetx,offsety,pval,px,py,stddev)
				}
				
				DLGvalue(called_dialog.lookupelement("driftx"),sum(out[0,0,1,sz])/sz)	//updating the UI fields with the mean values
				DLGvalue(called_dialog.lookupelement("drifty"),sum(out[1,0,2,sz])/sz)
				DLGvalue(called_dialog.lookupelement("driftphase"),sum(out[2,0,3,sz])/sz)
				
				//since this part is very prone to bugs it only operates on the output image. for proper drift corr the drift/correct buttons are intended to be used.
				//this section creates a lot of output for little gain, so deactivating it for now. add it in if you want to see a coarse drift-correction of the output collection
				
				//image reference_slices = realimage("",4,((sz-1)/2),7)
				//reference_slices = warp(out,(icol*2)+1,irow)		//reference slices are 2n+1
				//image measurement_slices = realimage("",4,((sz-1)/2),7)
				//measurement_slices = warp(out,(icol*2),irow)		//measurement slices are 2n
				//measurement_slices[0,0,3,((sz-1)/2)] += reference_slices[0,0,3,((sz-1)/2)] 	//correcting driftx/y/phase of measurement vs. immediate reference by reference vs. last reference
																								//+= instead of -= since values are flipped; should be verified
				//number avrgdriftx = sum(reference_slices[0,0,1,((sz-1)/2)])/((sz-1)/2)		//average values; using the pair-to-pair drift might be more precise
				//number avrgdrifty = sum(reference_slices[1,0,2,((sz-1)/2)])/((sz-1)/2)
				//number avrgdriftphase = sum(reference_slices[2,0,3,((sz-1)/2)])/((sz-1)/2)
				//DLGvalue(called_dialog.lookupelement("driftx"),avrgdriftx)
				//DLGvalue(called_dialog.lookupelement("drifty"),avrgdrifty)
				//DLGvalue(called_dialog.lookupelement("driftphase"),avrgdriftphase)
				
				//reference_slices.setname("0-tilt reference slices")	//the phase drift is very prone to 2pi() errors
				//showimage(reference_slices)
				//measurement_slices.setname("Corrected measurement slices")
				//showimage(measurement_slices)
			}
			
			profilestack.setname("Stack of extracted profiles")
			out.setname("Evaluation output")
			out.ImageSetDescriptionText("Xcorr X\nXcorr Y\nCoarse phase value\nSideband vector X\nSideband vector Y\nStd. dev.\nRefined phase value")
			showimage(out)
			
			if(showprofile) showimage(profilestack);									//the (real-space) profile stack extracted from the ROI for further processing
			
			image meas_profile_stack := realimage("Phase value / std.dev.",4,sz,1,2)	//we also choose to display the profile of the coarse phase + std.dev., since they are useful metrics for biprism drift and wobbler functionality 
			slice1(meas_profile_stack,0,0,0,0,sz,1) = slice1(out,0,2,0,0,sz,1)			//the final column (sz) has no meaning if double == 1, since it is the comparison of the last acquisition slice against itself
			slice1(meas_profile_stack,0,0,1,0,sz,1) = slice1(out,0,5,0,0,sz,1)
			showimage(meas_profile_stack)
			
			setpersistentnumbernote("Tilt Series:LastID", GetImageID(profilestack))		//store IDs in global tags to transfer images by reference later 
			setpersistentnumbernote("Tilt Series:OutID", GetImageID(out))
			TagGroup tgImg = ImageGetTagGroup(profilestack)
			if(!TagGroupDoesTagExist( tgImg, "PInit" ))	TagGroupCreateNewLabeledTag( tgImg, "PInit" )
			tgImg.TagGroupSetTagAsNumber( "PInit", SBLength )							//the mean 1*pi() wavelength is added as a tag to profile stack, to be used as an initial guess for the advanced fitting routine	
		}
		
		catch{	result("\nEvaluation ran into an exception."); called_dialog.halt_acquisition(); }	
		called_dialog.halt_acquisition()
		
		}
	
		if(switch == 2){	//case 2 is refine Xcorr measurement with the PCF method
							//the PCF currently does *not* respect the ROI 
		try{
		
			image current_out, Xcorr
			number outID
			getpersistentnumbernote("Tilt Series:OutID", outID)
			if(!GetImageFromID(current_out,outID)){ result("\nLast evaluation output not found!"); called_dialog.halt_acquisition(); exit(0); }
			
			get3dsize(front,sx,sy,sz)
			first = slice2(front,0,0,0,0,sx,1,1,sy,1)								//use the 0 slice as fixed reference, could also be the last slice (tiltback)
			
			self.Find_sideband(realFFT(first),px,py) 								//sideband pixel position in FFT space
			image aperture = self.Bworth_aperture(first,sqrt(px**2+py**2)*0.33,5) 	//hardcoded 5th order edge, 0.33x sideband distance, adjusting to a sharper edge/larger radius may result in better drift correction for some cases 
			//image aperture = self.Bworth_aperture(first[t,l,b,r],sqrt(px**2+py**2)*0.33,5)	//use this instead if you want to have ROI functionality in the PCF
			
			for(number i = 0; i<sz; i++){											//loop over every slice incl 0 (autocorrelation) to overwrite the previous results
				OpenAndSetProgressWindow( "Working on:", i+" of ", ""+sz )
				current = slice2(front,0,0,i,0,sx,1,1,sy,1)
				try{ Xcorr = self.Xcorr_norm(first,current,px,py,aperture); }		//measure drift vs first, could be done pairwise as well. not too hard to implement.
				//try{ Xcorr = self.Xcorr_norm(first[t,l,b,r],current[t,l,b,r],px,py,aperture); }	//use this instead if you want to have ROI functionality in the PCF
				catch{ px = 0; py = 0; }
				current_out.setpixel(i,0,(-1)*px)	//check the sign
				current_out.setpixel(i,1,(-1)*py)
			}
			current_out.setname("Evaluation output with PCF")
			
			DLGvalue(called_dialog.lookupelement("driftx"),sum(current_out[0,0,1,sz])/sz)	//updating the UI fields with the mean values
			DLGvalue(called_dialog.lookupelement("drifty"),sum(current_out[1,0,2,sz])/sz)
			DLGvalue(called_dialog.lookupelement("driftphase"),sum(current_out[2,0,3,sz])/sz)
		}
		
		catch{	result("\nDrift measurement ran into an exception."); called_dialog.halt_acquisition(); }	
		called_dialog.halt_acquisition()
		}
		
		if(switch == 3){	//case 3 is apply the current drift values to the stack; better to use DriftCorr_UI for detailed control
		try{
		
			image current_out
			number outID
			getpersistentnumbernote("Tilt Series:OutID", outID)
			if(!GetImageFromID(current_out,outID)){ result("\nLast evaluation output not found!"); called_dialog.halt_acquisition(); exit(0); }
			
			get3dsize(front,sx,sy,sz)
			image drift_out := imageclone(front)		//cloning the image stack can quickly lead to memory issues

			for(number i = 0; i<sz; i++){				//also align the first image vs itself, shouldnt make a difference since the autocorrelation maximum is always at 0,0
				current = slice2(front,0,0,i,0,sx,1,1,sy,1)
				OpenAndSetProgressWindow( "Working on:", i+" of ", ""+sz )
				px = getpixel(current_out,i,0)
				py = getpixel(current_out,i,1)
				try{ slice2(drift_out,0,0,i,0,sx,1,1,sy,1) = self.FFT_shift(current,px,py); }	//check the sign
				catch{ slice2(drift_out,0,0,i,0,sx,1,1,sy,1) = current; }
			}
			
			number maxsx = max(abs(current_out[0,0,1,sz]))	//unsigned max of drift components
			number maxsy = max(abs(current_out[1,0,2,sz]))
			
			image out_cut = drift_out[ceil(abs(maxsy)),ceil(abs(maxsx)),sy-ceil(abs(maxsy)),sx-ceil(abs(maxsx))]	//we cut the edges maximally, sensibly we would need to use the signed values t = max(drifty), l = max(driftx), b = sx-min(drifty), r = sy-min(driftx), I think (?)
			
			out_cut.ImageCopyCalibrationFrom(front)
			TagGroupCopyTagsFrom(imagegettaggroup(out_cut),imagegettaggroup(front))
			out_cut.setname(getname(front)+"_aligned_cut")
			showimage(out_cut)
		}
		
		catch{	result("\nApplying drift correction ran into an exception."); called_dialog.halt_acquisition(); }	
		called_dialog.halt_acquisition()
		
		}
	}
	
}

//constructor/destructor
evaluation_thread(object self)
{
}

~evaluation_thread(object self)
{
}
	
}

//building the UI
class Dialog_UI: UIFrame
{
object acquisition_thread, evaluation_thread //eventually we could also add a thread linking to Holoworks, for now this script is focused only on the acquisition of PS-EH stacks

number startangle,stopangle,stepangle	//UI field variables
number startamp,stopamp,stepamp
number driftx,drifty,driftphase
number xcal, ycal

//get the current user image shift via EM function, maybe also acquiring a reference image would be sensible? work in progress.
void toggleref(object self)
{
	number shiftx, shifty
	EMGetImageShift( shiftX, shiftY )
	result ("\nReference image shift x/y is : " + shiftx + "/" + shifty + "\n")
	setpersistentnumbernote("Tilt Series:ShiftX", shiftX)
	setpersistentnumbernote("Tilt Series:ShiftY", shiftY)
	self.halt_acquisition()
}

//initialize the threads and link them to the UI
object init(object self, number threadID)
{
	acquisition_thread = alloc(acquisition_thread)
	acquisition_thread.link_thread(threadID)
	evaluation_thread = alloc(evaluation_thread)
	evaluation_thread.link_thread(threadID)
	return self
}

//the UI reset
void halt_acquisition(object self)
{
	self.SetElementIsEnabled("startbutton",1)
	self.SetElementIsEnabled("wobblebutton",1)
	self.SetElementIsEnabled("stopbutton",0)
	self.SetElementIsEnabled("measurebutton",1)
	self.SetElementIsEnabled("refinebutton",1)
	self.SetElementIsEnabled("driftbutton",1)
	self.SetElementIsEnabled("correctbutton",1)
	do_break = 1
	exit(0)
	//no elegant way to halt thread execution except with the break flag or with setting signals.. currently threads simply escape the while() loop
}

//get field values and launch the acquisition thread to record PS-EH image stacks
//for a detailed description of the proper field values, also see the online documentation
void toggleon(object self)
{
	//get the current field values
	number pairwise = DLGGetValue(self.lookupelement("double"))
	number guntilt = DLGGetValue(self.lookupelement("guntilt"))
	number reference = DLGGetValue(self.lookupelement("refswitch"))
	
	startangle = DLGGetValue(self.lookupelement("anglestart"))
	stopangle = DLGGetValue(self.lookupelement("anglestop"))
	stepangle = DLGGetValue(self.lookupelement("anglestep"))
	
	startamp = DLGGetValue(self.lookupelement("ampstart"))
	stopamp = DLGGetValue(self.lookupelement("ampstop"))
	stepamp = DLGGetValue(self.lookupelement("ampstep"))
	
	//the angle values can be allowed to wrap around 2*pi, angle to run clockwise or counterclockwise, i think. 
	//still, cleaner to do it like this for now, means we cannot scan through a sector that e.g. crosses from -350 to +10 
	if((startangle>stopangle && stepangle>0) || stopangle<0 || (startamp>stopamp && stepamp>0) || stopamp <0 || stepangle<0 || stepamp<0){ OKDialog("\nReview your input values, currently only positive values and directions are allowed. \nWhen recording a static series, set both steps to 0 and fill the stop amp field with the number of iterations"); do_break = 1; return; }
	if(stepangle != 0 && stepamp != 0){ OKDialog("Only ramp either amplitude or angle."); do_break = 1; return; }
	
	acquisition_thread.init(pairwise,guntilt,reference,xcal,ycal,startangle,stopangle,stepangle,startamp,stopamp,stepamp)	//init the parameters
	
	//block the UI buttons
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	self.SetElementIsEnabled("stopbutton",1)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	
	do_break =0
	acquisition_thread.StartThread()
	//the choice of acquisition sequence is then made from within the acquisition thread depending on the switches and step fields set
}

//"wobble" button response, get the field values and launch the acquisition thread to record a wobbler series
void wobble(object self)
{
	number guntilt = DLGGetValue(self.lookupelement("guntilt"))
	if(!guntilt){ OKDialog("Only wobbling gun tilt supported!"); do_break = 1; self.halt_acquisition(); } 	//we currently only use the gun tilt API, the EMSetBeamTilt function is too slow for fast wobbling
																											//in order to wobble the beam tilt, you will need to write and compile your own API function
	number pairwise = 3												//we use the pairwise flag to call the wobbler functionality, since "pairwise wobbling" doesn't exist
	
	startangle = DLGGetValue(self.lookupelement("anglestart"))		//fixes the wobbler direction
	stopangle = DLGGetValue(self.lookupelement("anglestop"))		//is ignored
	stepangle = DLGGetValue(self.lookupelement("wobblesamples"))	//fixes the wobbling sampling points, passes into stepangle (dirty coding)
	
	startamp = DLGGetValue(self.lookupelement("ampstart"))			//fixes the wobbler minimum amplitude
	stopamp = DLGGetValue(self.lookupelement("ampstop"))			//fixes the wobbler maximum amplitude
	stepamp = DLGGetValue(self.lookupelement("ampstep"))			//fixes the number of steps between
	
	if(stopamp<startamp){ stopamp = startamp; startamp = 0; } //assuming user input error we wobble from 0 to max(start,stop) instead of trying to catch negative ampstep directions
	
	acquisition_thread.init(pairwise,guntilt,0,xcal,ycal,startangle,stopangle,stepangle,startamp,stopamp,stepamp)	//refswitch is static 0 since it doesn't have an influence on the wobbler
	
	//block the UI buttons
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	self.SetElementIsEnabled("stopbutton",1)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	
	do_break =0
	acquisition_thread.StartThread()
}

//"stop" button response, attempts to interrupt the threads
void toggleoff(object self)
{
	self.halt_acquisition()
}

//manually input calibration values, if you wish to use precise calibration values you will need to measure them via diffraction shift and store them in the global tags
void togglecalibrate(object self)
{
	xcal = GetNumber("X-tilt calibration value :",xcal,xcal)
	ycal = GetNumber("Y-tilt calibration value :",ycal,ycal)
}

//reworked drift measurement from center band PCF
void toggledrift(object self)
{
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	
	evaluation_thread.init(getfrontimage(),2,driftx,drifty,driftphase) //eval thread flag1 = 2 triggers reworked PCF drift measurement
																		//since the init() uses getfrontimage(), you will need to select the image stack again
	evaluation_thread.StartThread()
}

//quick measure button response, for coarse estimates of xcorr drift / phase shift etc.
void togglemeasure(object self)
{
	driftx = DLGGetValue(self.lookupelement("driftx"))
	drifty = DLGGetValue(self.lookupelement("drifty"))
	driftphase = DLGGetValue(self.lookupelement("driftphase"))
	
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	
	number pairwise = DLGGetValue(self.lookupelement("double"))
	
	evaluation_thread.init(getfrontimage(),pairwise,driftx,drifty,driftphase)
	do_break = 0
	evaluation_thread.StartThread()
	//run the evaluation thread over an entire stack, with flag1 = 0 for single and = 1 for pairwise with coarse (xcorr) drift correction
	//again the front image is assumed to be the acquired stack
}

//drifft correction of the current front image stack, using the values stored in the current OutID
void togglecorrect(object self)
{
	driftx = DLGGetValue(self.lookupelement("driftx"))		//the average drift values are outdated, leaving it here for possible later use
	drifty = DLGGetValue(self.lookupelement("drifty"))
	driftphase = DLGGetValue(self.lookupelement("driftphase"))
	
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	
	evaluation_thread.init(getfrontimage(),3,driftx,drifty,driftphase)
	do_break = 0
	evaluation_thread.StartThread()
}

//calls a script file in DMScriptPath, which in turn grabs the global tag values and runs a Python script
void togglerefine(object self)
{
	driftx = DLGGetValue(self.lookupelement("driftx"))
	drifty = DLGGetValue(self.lookupelement("drifty"))
	driftphase = DLGGetValue(self.lookupelement("driftphase"))
	
	self.SetElementIsEnabled("startbutton",0)
	self.SetElementIsEnabled("wobblebutton",0)
	
	self.SetElementIsEnabled("measurebutton",0)
	self.SetElementIsEnabled("refinebutton",0)
	self.SetElementIsEnabled("driftbutton",0)
	self.SetElementIsEnabled("correctbutton",0)
	//number pairwise = DLGGetValue(self.lookupelement("double"))	//currently not actually passed through
	
	//Version check for DM > 3.40 with Python functionality; of course, this depends on you having installed the Miniconda3 environment for DM... 
	number major, minor, debug
	getApplicationVersion(major,minor,debug)
	if(!((major > 2) && (minor >= 40))){ result("\nDM version too low for Python fitting."); self.halt_acquisition(); }
	else{ ExecuteScriptFile( DMScriptPath ); }		//see the CallPythonScript_mod.s file for suggestions on how to proceed from here
	//ideally we would like to try and write the results into the refined pval row, which means we need to find the output ID from the Python script... these are beta features that are up to the user to develop further
	self.halt_acquisition()
}

//build the UI
taggroup CreateDialog_UI(object self)
{
	startangle = 90; stopangle = 0; stepangle = 0;	//some default parameters in degree / [DAC] values
	startamp = 0; stopamp = 25; stepamp = 1;
	driftx = 0; drifty = 0; driftphase = 0;
	number samples = 2;
	xcal = 1; ycal = 1; //adjust once you have proper calibrated values, would need to distinguish between gun and beam

	taggroup dialog_items, label
	taggroup Dialog_UI = DLGCreateDialog("Tilt series UI", dialog_items)
	
	taggroup acquire_items
	taggroup acquire_box = DLGCreateBox("Acquire",acquire_items)
	taggroup startbutton = DLGcreatePushButton("Start","toggleon").DLGIdentifier("startbutton")
	taggroup refbutton = DLGcreatePushButton("Ref.","toggleref").DLGIdentifier("refbutton").DLGenabled(1)
	taggroup stopbutton = DLGcreatePushButton("Stop","toggleoff").DLGIdentifier("stopbutton").DLGenabled(0)
	taggroup topgroup = DLGGroupItems(startbutton, stopbutton, refbutton).dlgtablelayout(3,1,0)
	taggroup wobblebutton = DLGcreatePushButton("Wobble","wobble").DLGIdentifier("wobblebutton")
	label = DLGCreateLabel("Samples")
	taggroup wobblefield = DLGCreateRealField(samples,8,4).DLGIdentifier("wobblesamples")
	taggroup wobblegroup = DLGGroupItems(label, wobblefield).dlgtablelayout(1,2,0)
	taggroup wobbleunit = DLGGroupItems(wobblebutton,wobblegroup).dlgtablelayout(2,1,0)
	
	taggroup acquire_group = DLGGroupItems(topgroup,wobbleunit).dlgtablelayout(1,2,0).DLGinternalpadding(15,2)
	acquire_items.DLGAddElement(acquire_group)
	Dialog_UI.DLGAddElement(acquire_box)

	taggroup angle_items
	taggroup angle_box = DLGCreateBox("Angle [deg]",angle_items)
	label = DLGCreateLabel("Start")
	taggroup startangle_field = DLGCreateRealField(startangle,8,4).DLGIdentifier("anglestart")
	taggroup startangle_group = DLGGroupItems(label, startangle_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Stop")
	taggroup stopangle_field = DLGCreateRealField(stopangle,8,4).DLGIdentifier("anglestop")
	taggroup stopangle_group = DLGGroupItems(label, stopangle_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Stepsize")
	taggroup stepangle_field = DLGCreateRealField(stepangle,8,4).DLGIdentifier("anglestep")
	taggroup stepangle_group = DLGGroupItems(label, stepangle_field).dlgtablelayout(1,2,0)
	
	taggroup angle_group = DLGGroupItems(startangle_group,stopangle_group,stepangle_group).dlgtablelayout(3,1,0)
	angle_items.DLGAddElement(angle_group)
	Dialog_UI.DLGAddElement(angle_box)
	
	taggroup amp_items
	taggroup amp_box = DLGCreateBox("Amplitude [DAC]",amp_items)
	label = DLGCreateLabel("Start")
	taggroup startamp_field = DLGCreateRealField(startamp,8,4).DLGIdentifier("ampstart")
	taggroup startamp_group = DLGGroupItems(label, startamp_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Stop/Its.")
	taggroup stopamp_field = DLGCreateRealField(stopamp,8,4).DLGIdentifier("ampstop")
	taggroup stopamp_group = DLGGroupItems(label, stopamp_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Stepsize")
	taggroup stepamp_field = DLGCreateRealField(stepamp,8,4).DLGIdentifier("ampstep")
	taggroup stepamp_group = DLGGroupItems(label, stepamp_field).dlgtablelayout(1,2,0)
	
	taggroup amp_group = DLGGroupItems(startamp_group,stopamp_group,stepamp_group).dlgtablelayout(3,1,0)
	amp_items.DLGAddElement(amp_group)
	Dialog_UI.DLGAddElement(amp_box)
	
	taggroup options_items
	taggroup options_box = DLGCreateBox("Options",options_items)
	taggroup calibrate_button = DLGcreatePushButton("Calibrate","togglecalibrate").DLGIdentifier("calibratebutton").DLGenabled(1)
	taggroup gun_switch = DLGCreateCheckBox("Gun",1).DLGIdentifier("guntilt").DLGenabled(1)
	taggroup double_switch = DLGCreateCheckBox("Pairs",0).DLGIdentifier("double").DLGenabled(1)
	taggroup options_group = DLGGroupItems(calibrate_button,gun_switch,double_switch).dlgtablelayout(3,1,0)
	taggroup ref_switch = DLGCreateCheckBox("Vacuum reference",0).DLGIdentifier("refswitch").DLGenabled(1)	//switch to grab vacum references at the saved user image shift position
	
	options_items.DLGAddElement(options_group)
	options_items.DLGAddElement(ref_switch)
	Dialog_UI.DLGAddElement(options_box)
	
	taggroup eval_items
	taggroup eval_box = DLGCreateBox("Evaluate",eval_items)
	taggroup measurebutton = DLGcreatePushButton("Measure stack","togglemeasure").DLGIdentifier("measurebutton")
	taggroup driftbutton = DLGcreatePushButton("Drift","toggledrift").DLGIdentifier("driftbutton")
	taggroup eval_group = DLGGroupItems(measurebutton, driftbutton).dlgtablelayout(2,1,0).DLGinternalpadding(12,2)
	
	taggroup correctbutton = DLGcreatePushButton("Correct","togglecorrect").DLGIdentifier("correctbutton")
	taggroup refinebutton = DLGcreatePushButton("Refine","togglerefine").DLGIdentifier("refinebutton")
	taggroup eval_group2 = DLGGroupItems(correctbutton,refinebutton).dlgtablelayout(2,1,0).DLGinternalpadding(12,2)
	
	taggroup drift_items
	label = DLGCreateLabel("Drift X")
	taggroup driftx_field = DLGCreateRealField(driftx,8,4).DLGIdentifier("driftx")
	taggroup driftx_group = DLGGroupItems(label, driftx_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Drift Y")
	taggroup drifty_field = DLGCreateRealField(drifty,8,4).DLGIdentifier("drifty")
	taggroup drifty_group = DLGGroupItems(label, drifty_field).dlgtablelayout(1,2,0)
	label = DLGCreateLabel("Phase")
	taggroup driftphase_field = DLGCreateRealField(driftphase,8,4).DLGIdentifier("driftphase")
	taggroup driftphase_group = DLGGroupItems(label, driftphase_field).dlgtablelayout(1,2,0)
	taggroup drift_group = DLGGroupItems(driftx_group,drifty_group,driftphase_group).dlgtablelayout(3,1,0)
	
	eval_items.DLGAddElement(eval_group)
	eval_items.DLGAddElement(eval_group2)
	eval_items.DLGAddElement(drift_group)
	Dialog_UI.DLGAddElement(eval_box)
	
	taggroup footer=dlgcreatelabel(version)
	Dialog_UI.DLGAddElement(footer)
	
	taggroup position			
	position = DLGBuildPositionFromApplication()
	position.TagGroupSetTagAsString( "Width", "Medium" )
	position.TagGroupSetTagAsString( "Side", "Top" )
	Dialog_UI.DLGPosition(position) 

	return Dialog_UI
}

//constructor / destructor
Dialog_UI(object self)
{
	self.super.init(self.CreateDialog_UI())
}

~Dialog_UI(object self)
{	
	//result("\nTilt dialog interface destroyed.")
}

}

//main builds the UI dialog and creates the links
void main()
{
	object dialog_UI = alloc(Dialog_UI)
	Dialog_UI.display("Tilt series UI")
	Dialog_UI.init(Dialog_UI.ScriptObjectGetID())
}

main()
exit(0)