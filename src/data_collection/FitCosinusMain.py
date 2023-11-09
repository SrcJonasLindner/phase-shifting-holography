# Written by Jonas Lindner
# Solver Class adapted from: 
#https://omyllymaki.medium.com/gauss-newton-algorithm-implementation-from-scratch-55ebe56aac2e 


# Solver Class imports: 
from typing import Callable
import numpy as np
from numpy.linalg import pinv

# fit_data imports:
from numpy import exp, sin, cos, loadtxt
from matplotlib import pyplot as plt

# main imports: 
import DigitalMicrograph as DM

# Solver Class def ######################################################

class GNSolver:
    """
    Gauss-Newton solver.
    Given response vector y, dependent variable x and fit function f, 
    Minimize sum(residual^2) where residual = f(x, coefficients) - y.
    """

    def __init__(self,
                 fit_function: Callable,
                 max_iter: int = 1000,
                 tolerance_difference: float = 10 ** (-16),
                 tolerance: float = 10 ** (-9),
                 init_guess: np.ndarray = None,
                 ):
        """
        :param fit_function: Function that needs be fitted; y_estimate = fit_function(x, coefficients).
        :param max_iter: Maximum number of iterations for optimization.
        :param tolerance_difference: Terminate iteration if RMSE difference between iterations smaller than tolerance.
        :param tolerance: Terminate iteration if RMSE is smaller than tolerance.
        :param init_guess: Initial guess for coefficients.
        """
        self.fit_function = fit_function
        self.max_iter = max_iter
        self.tolerance_difference = tolerance_difference
        self.tolerance = tolerance
        self.coefficients = None
        self.x = None
        self.y = None
        self.init_guess = None
        if init_guess is not None:
            self.init_guess = init_guess

    def fit(self,
            x: np.ndarray,
            y: np.ndarray,
            init_guess: np.ndarray = None) -> np.ndarray:
        """
        Fit coefficients by minimizing RMSE.
        :param x: Independent variable.
        :param y: Response vector.
        :param init_guess: Initial guess for coefficients.
        :return: Fitted coefficients.
        """

        self.x = x
        self.y = y
        if init_guess is not None:
            self.init_guess = init_guess

        if init_guess is None:
            raise Exception("Initial guess needs to be provided")

        self.coefficients = self.init_guess
        rmse_prev = np.inf
        for k in range(self.max_iter):
           
            residual = self.get_residual()
            jacobian = self._calculate_jacobian(self.coefficients, step=10 ** (-6))
            #print('calc jac: ' , jacobian[0:3,0:5]) 
            
            self.coefficients = self.coefficients - self._calculate_pseudoinverse(jacobian) @ residual
     
            #print('calc coeff: ' , (self._calculate_pseudoinverse(jacobian) @ residual)) 
            #print('shape pseudo: ' , np.shape(self._calculate_pseudoinverse(jacobian)) )
            #print('calc pseudo: ' , (self._calculate_pseudoinverse(jacobian)[0:5,0:3])) 
             
            #print("pseudo : ", self._calculate_pseudoinverse(jacobian)[0:3,52], "\n")
            rmse = np.sqrt(np.sum(residual ** 2))
            #logger.info(f"Round {k}: RMSE {rmse}")
            if self.tolerance_difference is not None:
                diff = np.abs(rmse_prev - rmse)
                if diff < self.tolerance_difference:
                    #logger.info("RMSE difference between iterations smaller than tolerance. Fit terminated.")
                    #RMSEs.append(rmse)
                   #print('RMSE: ', rmse)
                    return self.coefficients
            if rmse < self.tolerance:
                #logger.info("RMSE error smaller than tolerance. Fit terminated.")
                #RMSEs.append(rmse)
                #print('RMSE: ', rmse)
                return self.coefficients
            rmse_prev = rmse
        #print("Max number of iterations reached. Fit didn't converge!")   
    
        return np.array([-1,-1,-1,-1])

    def predict(self, x: np.ndarray):
        """
        Predict response for given x based on fitted coefficients.
        :param x: Independent variable.
        :return: Response vector.
        """
        return self.fit_function(x, self.coefficients)

    def get_residual(self) -> np.ndarray:
        """
        Get residual after fit.
        :return: Residual (y_fitted - y).
        """
        return self._calculate_residual(self.coefficients)

    def get_estimate(self) -> np.ndarray:
        """
        Get estimated response vector based on fit.
        :return: Response vector
        """
        return self.fit_function(self.x, self.coefficients)

    def _calculate_residual(self, coefficients: np.ndarray) -> np.ndarray:
        y_fit = self.fit_function(self.x, coefficients)
        return y_fit - self.y

    def _calculate_jacobian(self,
                            x0: np.ndarray,
                            step: float = 10 ** (-6)) -> np.ndarray:
        """
        Calculate Jacobian matrix numerically.
        J_ij = d(r_i)/d(x_j)
        """
        y0 = self._calculate_residual(x0)

        jacobian = []
        for i, parameter in enumerate(x0):
            x = x0.copy()
            x[i] += step
            #if i==0:
                #print(x)
            
            y = self._calculate_residual(x)
            derivative = (y - y0) / step
            jacobian.append(derivative)
        jacobian = np.array(jacobian).T

        #print(jacobian[52,2])
        return jacobian

    @staticmethod
    def _calculate_pseudoinverse(x: np.ndarray) -> np.ndarray:
        """
        Moore-Penrose inverse.
        """
        return pinv(x.T @ x) @ x.T





# fit_data ############################################

def func(x, coeff):
    return coeff[0] * np.cos(coeff[1]*(x)+coeff[2])    
    
def fit_data(raw_data, T,plot_fits,phase_init_offset):
    
    [_,dim]=np.shape(raw_data)
    x= np.array(range(0,dim))
    data=raw_data[0][:]

    # remove offset and normalize data:
    data=data-np.mean(data)
    data = data/np.linalg.norm(data)

    start_amp= (max(data)-min(data))/2

    # 2*Pi()/ (T) = start freq
    #T=17.0 # period in pixel
    
    start_freq=(2.0*np.pi)/T
    start_phase=((np.argmax(data[1:int(T*1.5)], axis=0))/T*2.0*np.pi)%(2.0*np.pi)
    #print('max pos: ',np.argmax(data[1:int(T*1.5)], axis=0) )
    #print('start phase: ' , start_phase)
    #print((2*np.pi)/T)
    # tolerance diff 1e-9!
    
   
    COEFFICIENTS = [start_amp,start_freq ,start_phase+phase_init_offset]
    #print('start coeff: ' , COEFFICIENTS)
    
    solver = GNSolver(fit_function=func, max_iter=1000, tolerance_difference=10 ** (-9))
    init_guess = COEFFICIENTS
    _ = solver.fit(x, data, init_guess)
    fit = solver.get_estimate()
    residual = solver.get_residual()
    rmse = np.sqrt(np.sum(residual ** 2))
    
    R_sq = 1 - sum((residual)**2)/sum((data-np.mean(data))**2)
    
    
    # p == number of fitt vars
    p=np.shape(COEFFICIENTS)[0]
    adj_R_sq= 1 - (1 - R_sq) * ((dim - 1)/(dim-p-1))

    if (plot_fits==1):
        plt.figure()
        plt.plot(x, data, label="Original, noiseless signal", linewidth=2)
        plt.plot(x, fit, label="Fit", linewidth=2)
        plt.plot(x, residual, label="Residual", linewidth=2)
        plt.title("Gauss-Newton: curve fitting example")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid()
        plt.legend()
        plt.show()
        #print(solver.coefficients)
    
    fit_amp=   solver.coefficients[0]
    fit_omega= solver.coefficients[1]
    fit_phase= solver.coefficients[2]

    return [fit_amp,fit_omega,fit_phase,adj_R_sq];
  
# the main script ############################################

def main():

	print('starting fitting\n')
	
	#<InputStart>
	InputID=973
	print(InputID,' Input ID')
	T=17.7086
	plot_fits=0 # still buggy
	#<InputEnd>

	stack = DM.FindImageByID(InputID)
	stackTags = stack.GetTagGroup()
	StackData = stack.GetNumArray() 
	
	# Line Profile dims: [sz,sy,sx]
	[sz,sy,sx]= StackData.shape
	
	#print(StackData.shape, type(StackData))
	
	OutputData=[]
	
	for i in range(0,sz):
		print('fitting slice ', i+1 ,'/' , sz)
		#Paras:  [amp, freq, phase, adj. R^2]
		paras=[0.0,0.0,0.0,0.0]
		raw_data=np.array(StackData[i,:,:])
		
		phase_init_offset=0.0
		paras=fit_data(raw_data,T,plot_fits,phase_init_offset);
		
		
		k=0
		while (paras[3]<=0.75 and k<=16):
			#print('fit was poor, adjunsting phase init')
			phase_init_offset+=1.0/8.0*np.pi
			paras=fit_data(raw_data,T,plot_fits,phase_init_offset );
			k+=1
		
		if (paras[3]<0.75): print('Warning: Fit did not Converge on Slice ', i  ,"!\n")
		
		# if fitting amp negative: add Pi()
		if (paras[0]<=0): 
			paras[0]*=-1.0
			paras[2]+=np.pi
		paras[2]=paras[2]%(2.0*np.pi)
		OutputData.append(paras)
		
	OutputData=np.array(OutputData)
	OutputImg= DM.CreateImage(OutputData)

	
	OutputImg.ShowImage()
	ReturnID= OutputImg.GetID()
	OutputImgTags=OutputImg.GetTagGroup()
	
	OutputImgTags.SetTagAsString( 'Info:Description', " [amp,omega,phase,adj. R^2] x [sz]\n" )
	
	stackTags.SetTagAsLong( 'ReturnImageID', ReturnID )
	del stack
	del OutputImg
	
main()