# written by Jonas Lindner, Univeristy of Goettingen,  2022
# accordong to  https://doi.org/10.4028/www.scientific.net/MSF.833.215


#from ru et al 1994
# calculate and invert the matrix, that does a linear Fourier series fit for each pixel.
# Matrix inversion is done in python, since DM supports no double prec. matrix inversion at this time.
# This is a python scritpt should be started from Digital Micrograph. Select lenguage python in the bottom of this window.
# And of course instalkl pyhton support.

import DigitalMicrograph as DM
import numpy as np
from numpy import exp, complex

"""
Input: Takes a profile of phase extracted and analysed by TiltSeries_UI.s Measure button.
       The phase is located in the third row of the output image. Draw a line profile is necessary.
       The image must be opened in Digital mircograph and the Label passed as input.

Output: complex 3x3 inverted Matrix as two images: one contains the real and the other the imaginary part 
further processed with REcon PS holo image 

Check if input lineprofile contains empty pixels at start or end!
"""


"""
Input:  -Takes a profile of phase extracted and analysed by cosFitMain.py
		-Takes a vector of k values to include in the fourier fitting series. 
		 1 = sideband k; Keep the sceme if adding new values e.g. [0,-1.0,1.0,-0.5,0.5]! 
Output: complex nxn inverted Matrix as two images: one contains the real and the other the imaginary part 
further processed with REcon PS holo image 

Check if input lineprofile conains empty pixels at start or end!
"""

def main():
	#<InputStart>
	ProfileLabel="CW"
	vec_k=[0,-1.0,1.0,-0.5,0.5,-0.25,0.25] 
	#<InputEnd>
	
	order_L=len(vec_k)
	
	try: 
		profile= DM.FindImageByLabel(ProfileLabel)
		profiledata=profile.GetNumArray()
	except:
		print('Input Images not found!')
		exit(0)
		
	# para matrix: prefactor of Phi_n in eqation (8) in ref
	print(order_L)
	para_matrix=np.zeros((order_L,order_L))
		
	for m in range(0,order_L):
		for j in range(0,order_L):
			# swapped [m,j] to [j,m] due to python-> DM image convention
			para_matrix[j,m]= vec_k[j]-vec_k[m] # swapped j,m due to python-> DM image convention
	print(para_matrix)

	# building the matrix in eq (8)
	
	matrix=np.zeros((order_L,order_L), dtype=complex)
	dim = int(np.shape(profiledata)[0])
	print("dim=",dim)
	for i in range(0,order_L):
		for j in range(0,order_L):
			
			if (i==j): # diag element= N
				matrix[i,i]= dim;
			else: # off dia == sum(exp(paras_matrix[i,j] * i *phi0(n)) 
				val=complex(0.0,0.0)
				for k in range(0,dim):
					val+= np.exp(para_matrix[i,j]*complex(0.0,1.0)* profiledata[k]); 
				matrix[j,i]=val;
			
	invert= np.linalg.inv(matrix)
	print('inverse matrix: ')
	print(invert)

	invert_real= np.zeros((order_L,order_L))
	invert_real+= np.real(invert)
	
	invert_img= np.zeros((order_L,order_L))
	invert_img+=  np.imag(invert)
	
	
	
	imag= DM.CreateImage(invert_img)
	real= DM.CreateImage(invert_real)
	vec_k_img=DM.CreateImage(np.array(vec_k))
	
	real.SetName("real inverted matrix")
	imag.SetName("imag inverted matrix")
	vec_k_img.SetName("norm k values")
	
	vec_k_img.ShowImage()
	imag.ShowImage()
	real.ShowImage()
	
	del imag
	del real
	del vec_k_img
	del profile

main()