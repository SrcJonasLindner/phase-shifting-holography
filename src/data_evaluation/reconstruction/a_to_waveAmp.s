// normalize a and a_ref to Wave Amp A()
//  see Eq. (7) the publication. 

// labels of a and a_:ref
// note that the std output from Recon_Image.s is a complex image with only values in Re()
compleximage a:=AA
compleximage a_ref:=Y

image amp=real(a)/(real(a_ref)/2)-1

image counter=real(a) 
counter=0 // counter for problematic pixels, that have to be set to zero. Be alerted if >2.

counter =tert(amp<0 || isNan(amp),1,0)
amp= tert(amp<0 || isNan(amp),0,amp)

result("sum problematic pixels =" +sum(counter)+"\n")
amp.setname("sqrt(2a/a_ref-1)")
showimage(sqrt(amp))



