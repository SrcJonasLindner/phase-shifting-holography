// calculate R^2 and R^2 adj for a given PS-Reconstruction



image ps_stack= AQ
realimage fit_stack=real(J)

number sx,sy,sz, Tx, Ty, p=3.0 // p= number of parameters in fit

get3Dsize(ps_stack,sx,sy,sz)

image residual_sq_sum= realimage("",4,sx,sy)
image std_derv_data=realimage("",4,sx,sy)

image R_sq=  realimage("",4,sx,sy)
image R_sq_adj=  realimage("",4,sx,sy)

image mean_data= realimage("",4,sx,sy)

for (number i=0; i<sz; i++){
	image current_data= slice2(ps_stack, 0,0,i, 0, sx, 1, 1, sy,1)
	mean_data+=current_data/sz
}

// use the fit parameter to calc the fitted cos wave for each pixel projection
for (number i=0; i<sz;i++){

image current_data= slice2(ps_stack, 0,0,i, 0, sx, 1, 1, sy,1)
image current_fit= slice2(fit_stack, 0,0,i, 0, sx, 1, 1, sy,1)

// calc residual between fit and exp data:
residual_sq_sum+= (current_data-current_fit)**2

// calc stad abw of exp data:
std_derv_data+= (current_data-mean_data)**2

}

//R^2:
R_sq= 1.0-residual_sq_sum/std_derv_data
R_sq.setName("R_sq")
showimage(R_sq)
//adj. R^2:
R_sq_adj= 1.0 - (1.0 - R_sq) * ((sz - 1.0)/(sz-p-1.0))

R_sq_adj.setName("R_sq_adj")
showimage(R_sq_adj)

/* rmse = np.sqrt(np.sum(residual ** 2)) 
    R_sq = 1 - sum((residual)**2)/sum((data-np.mean(data))**2)
    p=3
    adj_R_sq= 1 - (1 - R_sq) * ((dim - 1)/(dim-p-1))
*/

