
// written by Jonas Lindner, Univeristy of Goettingen,  2022
// Phaseshifting reconstruction
// with fresnel correction by  Lei et al: 
//according to https://doi.org/10.4028/www.scientific.net/MSF.833.215

// implemets ru et al. 

// Tx/Ty ramp can be also removed from 'phase' image with post processing
// normalized k in vec_k= k(k_j-k_m*)/k_sideband
// vec(C_j) = matrix(M)^-1* vec(v_j)


number sx,sy,sz, Ty, Tx
image stack, real, imag, phaseFit, vec_k
try{
	 stack= C // stack with holograms aquired with phase shifting holo
	 phaseFit= F    // initial phase shifts of the holo-stack as profile
	 real= I        // real part of the inverted matrix eq. (8)
	 imag= B       // imag part of the inverted matrix eq. (8)
	 vec_k=J      // the normalized k_values used for Matrix recon as profile
	 Tx=36.0     // fringe period in x direction [px]
	 Ty=14.9     // fringe period in y direction [px]

	}
catch{
	result("missing at least one input image: stack, real mat, img mat, PhaseFit")
	exit(-1)
}
number order_L;
getsize(real,order_L,order_L)

get3DSize(stack, sx,sy,sz);


compleximage Cj_stack= compleximage("",8,sx,sy,order_L)
compleximage vj_stack= compleximage("",8,sx,sy,order_L)
compleximage matrix = complex(real,imag)


for (number j=0; j<order_L; j++){
	
	compleximage vj_act= compleximage("",8,sx,sy)
	for (number i=0; i<sz; i++){ // for each hologram in stack
		
		// vj implents right hand side of eq. 8 in Lei et al.
		number k_val= getPixel(vec_k,j,0) // normalized get k-value
		vj_act+= slice2(stack, 0,0,i, 0, sx, 1, 1, sy,1)*exp( k_val *complex(0.0,1.0)*getPixel(phaseFit,i,0))	
	} 
	slice2(vj_stack, 0,0,j, 0, sx, 1, 1, sy,1)=vj_act;
}

for (number j=0; j<order_L; j++){
	compleximage C_act= compleximage("",8,sx,sy)
	
	for (number k=0; k<order_L; k++){// C_j= sum(v_k * m_jk)
		
		C_act+= slice2(vj_stack, 0,0,k, 0, sx, 1, 1, sy,1)*getpixel(matrix,j,k);
	}
	slice2(Cj_stack, 0,0,j, 0, sx, 1, 1, sy,1)=C_act;
}

//still buggy
// calc the fit stack
compleximage fit_stack=compleximage("fit_stack",8,sx,sy,sz)
fit_stack.setName("fit_stack")

for (number i=0; i<sz; i++){
	// for each hologram
	compleximage fit_act= compleximage("",8,sx,sy)
	number act_phase=getpixel(phaseFit,i,0)
	// sum(exp(i*k_j)*C_j)
	for (number j=0; j<order_L; j++){
		number k_val=getPixel(vec_k,j,0);
		compleximage C_act=slice2(Cj_stack, 0,0,j, 0, sx, 1, 1, sy,1)
		fit_act+= C_act*exp(-complex(0,1)*k_val*act_phase)
	}
	slice2(fit_stack, 0,0,i, 0, sx, 1, 1, sy,1)=fit_act;
}

Cj_stack.setname("Cj_stack")
Cj_stack.ImageSetDescriptionText("Tx="+Tx+", Ty="+Ty)
//showimage(Cj_stack)

showimage(fit_stack)

/*
phase= atan(Imaginary(C2)/real(c2));

phase.SetName("phase")

// sign convention likely different than in paper, therefore add the Tx/Ty phaseramp
phase2=atan(imaginary(C2)/real(C2))
//phase2 = phase2+mod(2*pi()*(icol)/(Tx),pi())+mod(2*pi()*(irow)/(Ty),Pi())
//phase2 = phase2+mod(2*pi()*(icol)/(Tx)+2*pi()*(irow)/(Ty),Pi())

phase2 = mod(phase2-2*pi()*(icol)/(Tx)-2*pi()*(irow)/(Ty),Pi())
phase2.ImageSetDescriptionText("Tx="+Tx+", Ty="+Ty)

number eps=1e-3
phase2=tert(abs(phase2-pi())>eps,phase2,phase2-pi())



*/