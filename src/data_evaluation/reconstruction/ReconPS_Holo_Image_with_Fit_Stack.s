
// second part of the phase shifting holography reconstruction
// from Ru et al. 1994 ultra microscopy eq (8) part 2: inverse matrix * img_vec

/* ##### Inputs: (image labels, numbers)
	 stack= AG // experimental holo stack, specimen drift corrected with Drift_corr_UI
	 phaseFit= AU // line profile that contains the phase of the carrier freq (after drift corr)
	 real= AV // real part of the 3x3  inverted matrix as calc by Recon matrix.py
	 imag= AT // imag part of the 3x3  inverted matrix as calc by Recon matrix.py
	 Tx=46.5     // fringe period in x direction [px] Use the same parameters for ref. reconstruction.
	 Ty=18.79    // fringe period in y direction [px] Use the same parameters for ref. reconstruction.

	//Note that the sign of Tx and Ty depends on the choice of the side band maximum to calculate the phase shift. And may have to be adapted. 
 	//By default (both positive) we assume, that the one that is used is in the upper left quadrant. 
 	//If it is located in the upper right quadrant, flip the sign of $T_x$.
  
   ### Outputs:
   
   image a= offset of the linear cosine fit result for each pixel 
   image b= amplitude of the the linear cosine fit result for each pixel
   image phase= Phase fit without substraction of Tx/Ty ramp. Can be used to refine Tx,Ty
   image phase2=Phase fit with substraction of Tx,Ty ramp. 
   image fit_stack= stack of holograms calculated from the fit parameters. 
                    The difference to the experimental series can be used to calculate the goodness of fit (R^2 or R^2 adj.). 
*/


number sx,sy,sz, Ty, Tx
image stack, real, imag, phaseFit

// ### Inputs: 
try{ // insert the image labels/numbers  here
	 stack= AG // holo stack
	 phaseFit= AU // phase line profile
	 real= AV  // real part matrix 3x3
	 imag= AT  // imag part matrix 3x3
	 Tx=46.5   // finge period x
	 Ty=18.79  // finge period y  
	}
catch{
	result("missing at least one input image: stack, real mat, img mat, PhaseFit")
	exit(-1)
}

get3DSize(stack, sx,sy,sz);

compleximage fit_stack= compleximage("",8,sx,sy,sz)
compleximage C1= compleximage("",8,sx,sy)
compleximage C2= compleximage("",8,sx,sy)
compleximage C3= compleximage("",8,sx,sy)

compleximage v1= compleximage("",8,sx,sy)
compleximage v2= compleximage("",8,sx,sy)
compleximage v3= compleximage("",8,sx,sy)

compleximage matrix = complex(real,imag)

for (number i=0; i<sz; i++){
	v1+= slice2(stack, 0,0,i, 0, sx, 1, 1, sy,1)
	v2+= slice2(stack, 0,0,i, 0, sx, 1, 1, sy,1)*exp( (-1.0) *complex(0.0,1.0)*getPixel(phaseFit,i,0))
	v3+= slice2(stack, 0,0,i, 0, sx, 1, 1, sy,1)*exp(   1.0  *complex(0.0,1.0)*getPixel(phaseFit,i,0))
} 

C1= v1*getpixel(matrix,0,0)+ v2 * getPixel(matrix,0,1) + v3* getPixel(matrix,0,2);
C2= v1*getpixel(matrix,1,0)+ v2 * getPixel(matrix,1,1) + v3* getPixel(matrix,1,2);
C3= v1*getpixel(matrix,2,0)+ v2 * getPixel(matrix,2,1) + v3* getPixel(matrix,2,2);


// calc the fit stack
for (number i=0; i<sz;i++){

number act_phase=getPixel(phaseFit,i,0)
compleximage current= C1+ C2*exp(1.0 *complex(0.0,1.0)*act_phase) + C3*exp(-1.0 *complex(0.0,1.0)*act_phase)
slice2(fit_stack, 0,0,i, 0, sx, 1, 1, sy,1)=real(current)
}

fit_stack.SetName("fit_stack")
showimage(fit_stack)

compleximage a,b
image phase, phase2

a=C1
b=2*sqrt(C2*C3);
phase= atan2(Imaginary(C2),real(C2));

a.SetName("a")
b.SetName("b")
phase.SetName("phase")

phase2=atan2(imaginary(C2),real(C2))
phase2 = mod(phase2-2*pi()*(icol)/(Tx)-2*pi()*(irow)/(Ty),2*Pi())

phase2.SetName("phase2")
showimage(a)
showimage(b)
showimage(phase)
showimage(phase2)
