import DigitalMicrograph as DM
import numpy as np
from numpy import exp, complex
#from ru et al 1994
# calculate and invert the matrix, that does a linear cosine fit for each pixel.
# Matrix inversion is done in python, since DM supports no double prec. matrix inversion at this time.
# This is a python scritpt should be started from Digital Micrograph. Select lenguage python in the bottom of this window.
# And of course instalkl pyhton support.

"""
Input: Takes a profile of phase extracted and analysed by TiltSeries_UI.s Measure button.
       The phase is located in the third row of the output image. Draw a line profile is necessary.
       The image must be opened in Digital mircograph and the Label passed as input.

Output: complex 3x3 inverted Matrix as two images: one contains the real and the other the imaginary part 
further processed with REcon PS holo image 

Check if input lineprofile contains empty pixels at start or end!
"""

def main():
	#<InputStart>
	ProfileLabel="E"
	#<InputEnd>
	
	paras=np.array([[  0.,  1.0, -1.0],
					[ -1.0,  0., -2.0],
					[  1.0,  2.0, 0.]])
	
	profile=None;
	try:
		profile= DM.FindImageByLabel(ProfileLabel)
	except:
		print('Input Images not found!')
		exit(0)
	profiledata=profile.GetNumArray()

	print(np.shape(profiledata))
	
	matrix=np.zeros((3,3), dtype=complex)
	dim = int(np.shape(profiledata)[0])
	print("dim=",dim)
	for i in range(0,3):
		for j in range(0,3):
			
			if (i==j): # diag element= N
				matrix[i,i]= dim;
			else: # off dia == sum(exp(paras[i,j] * i *phi0(n)) 
				val=complex(0.0,0.0)
				for k in range(0,dim):
					val+= np.exp(paras[i,j]*complex(0.0,1.0)* profiledata[k]); 
				matrix[j,i]=val;
			
	invert= np.linalg.inv(matrix)
	print('inverse matrix: ')
	print(invert)

	invert_real= np.zeros((3,3))
	invert_real+= np.real(invert)
	
	invert_img= np.zeros((3,3))
	invert_img+=  np.imag(invert)
	
	imag= DM.CreateImage(invert_img)
	real= DM.CreateImage(invert_real)
	
	imag.ShowImage()
	real.ShowImage()
	
	imag.SetName("imaginary invert matrix")
	real.SetName("real invert matrix")
	
	del profile
	del real
	del imag

main()