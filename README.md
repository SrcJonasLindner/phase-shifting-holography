# Phase-shifting-holography software package

[![Phase Shifting Paper](https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/PaperHeader.png)](https://doi.org/10.1016/j.ultramic.2022.113629)

Software package collecting scripts for the data collection and evaluation of phase shifting holography image stacks with Digital Micrograph Scripting and TitanScripting for electron microscopy.

--- 
**Table of contents:**

1. [Requirements](#requirements)
1. [Data Collection](#data-collection) <br>
    1. [Tilt Series UI](#tilt-series-ui) <br>
    1. [Calibrating the tilt](#calibrating-the-tilt)
1. [Data Evaluation](#data-evaluation)<br>
    1. [Specimen drift correction](#specimen-drift-correction)<br>
    1. [Phase shifting series reconstruction](#phase-shifting-series-reconstruction)<br>   
        1. [Measure the carrier frequency phase](#measure-the-carrier-frequency-phase)<br>
        1. [Reconstruction Matrix](#reconstruction-matrix)<br>
        1. [Reconstruction Image](#reconstruction-image)<br>
1. [Example Data and Workflow](#example-data-and-workflow)<br>
1. [Copyright and Diclaimer](#copyright-and-disclaimer)
---

# Requirements

  ## Digital Micrograph 3.5X+ with Python support installed
  
  Install the latest version [Digital Micrograph](https://www.gatan.com/installation-instructions) version
  
  ## Software Microscope control to lenses (optional, but strongly recommended)
  The phase shifting measurement procedure involes a routine to use the GunTilt lens of the Titan to tilt the beam off from the optical axis to create phase    shifts in the holograms, that are the heart of this method. There is a [default option](http://www.dmscripting.com/tem_control.html) using the Digital Micrograph Microscope interface, but with our Microscope this procedure was way too slow and needed approx 1-2 seconds to create one of the beam tilts. Since a phase shifting holo series covers around several tens of images, the specimen driftin the dead time is a serious problem. 
 
  [Tem scripting interface](https://temscript.readthedocs.io/) or the ability to precompile microscope commands as executable files, that then are included in the measurement script.
  
---
# Data collection

The acquisition of the phase shifting holo image series requires a calibration of the magnitude and direction of the beam tilt. Both, calibration and the image series recording are implemented in the [TiltSeriesUI.s]() Digital Micrograph script. An phase shifting holo tilt series should at least cover the phaseshifts in the intervall $\[0,2\pi\]$. Generally choose short exposure time (e.g. one second) to be able to correct the specimen drift via post processing. Aim for 50 holograms in the tilt series to increase the fitting accuracy of the reconstruction and high fringe visibility.
 
---
  ## Tilt Series UI

<div class="content">
  <img align="right" src="https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/TiltSeriesUI.png">
</div>

<details>
 <summary>Acquire</summary>
 
 <kbd>Start</kbd>: Starts the recording of the image tilt series. Looping over the amplitude and angle ranges given under<br>&emsp;&emsp;&emsp; [Amplitude](#tilt-series-ui) and [Angle](#tilt-series-ui).<br>
 <kbd>Stop</kbd>: Cancels the measurement thread.<br>
 <kbd>Ref</kbd>: Saves the current image shift position to jump there to take the vacuum reference after the specimen tilt<br>&emsp;&emsp; series.<br>
 <kbd>Wobble</kbd>: Starts a Wobbler-calibration measurement. Using  [Amplitude](#tilt-series-ui),  [Angle](#tilt-series-ui) and  [Samples](#tilt-series-ui) <br>
 <kbd>Samples</kbd>: Defines the number of sample-points within the wobbler amplitude ramp. See [details](#calibrating-the-tilt)

</details>



<details>
  <summary>Angle [deg]</summary>
    
 <kbd>Start</kbd>: Inital tilt angle. 90 degree corresponds to the tem holder axis in the FEI Titan.<br>
 <kbd>Stop</kbd>: Final tilt angle.<br>
 <kbd>Stepsize</kbd>: Used stepsize if ramping the angle. Set to zeroi for fixed angle measurements.<br>
  
</details>

<details>
  <summary>Amplitude [DAC]<a name="deg"></a></summary>
    
  <kbd>Start</kbd>: Initial tilt amplitude in units to the digital to analog converter.  <br>
 <kbd>Stop/Its.</kbd>: Final tilt magnitude. For typical biprism voltages from 200-300 V [0,40] is a good first wobbbling<br> &emsp;&emsp;&emsp;&emsp; guess.<br>
 <kbd>Stepsize</kbd>: If stepsize=0: using the stop value as number of measurements to be done. If stepsize>=0 use it to<br>&emsp;&emsp;&emsp;&emsp; determine the steps between <kbd>Start</kbd>and <kbd>Stop</kbd>.<br>
  
</details>

<details>
  <summary>Options</summary>
    
<kbd>Calibrate</kbd>: User input of constant factors to scale the tilt magnitudes in x- and y directions<br>&emsp;&emsp;&emsp;&emsp; respectively. Optional after calibration to rescale the used tilt magnitude.<br>    
[&check;] Gun: If checked the gun tilt lens is used via [TEM scripting](https://temscript.readthedocs.io/) compiled executable. <br>&emsp;&emsp;&emsp;&emsp;If unchecked the [EM Commands](http://www.dmscripting.com/tem_control.html) are used. Note that the later may have a long delay. <br>
[&check;] Pairs: If checked an image is recorded at the zero reference tilt angle after each image in the tilt series. This<br> &emsp;&emsp;&emsp;&emsp;feature can be used for [π-phase shifting electron holography](https://doi.org/10.1016/j.ultramic.2018.06.004) <br>
[&check;] Vacuum reference: If checked a vacuum reference tilt series is acquired parallel at the image shift position<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp; defined by <kbd>Ref.</kbd> The same angle and amplitude parameters are used for the vacuum<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp; series. Trades specimen drift vs. biprsim drift compared to sequenctial recording of two<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp; seperate stacks.<br>
    
</details>

<details>
  <summary>Evaluate</summary>
    
<kbd>Drift</kbd>: Uses the [phase correlation function](https://doi.org/10.1016/S0304-3991(02)00071-2) to roughly calculate the specimen drift  of the<br>&emsp;&emsp;&emsp;&nbsp; frontwindow <br>&emsp;&emsp;&emsp;&nbsp;image stack.<br>
<kbd>Measure stack</kbd>: Measures the hologram phase of the carrier frequency from the FFT of each hologram inside the frontmost stack. The result is provided as line profile, that is used as input for the [reconstruction](#phase-shifting-series-reconstruction).  <br>
<kbd>Correct</kbd>: Corrects the tilt series image stack by the drift determined by <kbd>Drift</kbd>.<br>
<kbd>Refine</kbd>: Beta feature that measures the phase of the holograms by non-linear cosine fitting via python.<br>&emsp;&emsp;&emsp;&nbsp; Not recommended, use <kbd>Measure stack</kbd> insted.<br>
    
<kbd>Drift x</kbd>: Linear x-drift interpolation according to <kbd>Drift</kbd>. <br>
<kbd>Drift y</kbd>: Linear y-drift interpolation according to <kbd>Drift</kbd>. <br>
</details>    


  ---
  ## Calibrating the tilt
    
For a given biprism voltage and a certain beam tilt angle the magnitude $\phi$ of the beam tilt in DAC-units has to be calibrated. To minimize the influence of mechanical vibration the biprism should be oriented perpendicular to the holder axis. In the Titan the coordinate system of the beam deflection the y-direction is parallel to the holder axis. The wobbling calibration procedure takes a series of images. The exposure time of each image is synchronized with a tilt-ramp during this exposure. The number of points within the ramp $n$ (samples) has to be chosen beforehand. Each image is acquired with a different maximum tilt $\phi_{max}$. The intensity of the acquired wobbling-image can be described as:
    
$$ I_{wobbler}(\phi_{max},n)=\sum_{i=1}^{n} a(x,y) + b(x,y) \cos \left[ \frac{2\pi x}{Tx} + \frac{2\pi y}{Ty} + \frac{\phi_{max}}{(n-i)}  \right]  $$

If the standard derivation of each image is plotted against the $\phi_{max}$ or tilt [DAC] the tilt magnitude needed to create a phase shift of $2\pi n$ can be measured by the distance of maxima of the standard derivative. The standard derivative is maximized if each step corresponds to $\phi_{max}/((n-i))=2\pi m$. This procedure allows a fast way of calibrating the needed tilt for arbitrary biprism orientations and angles and is included in the UI-script.
    
 ```
    Wobbler-Calibration Input Example (250V biprsim): 
    
    Acquire:
        Samples=4.0
    Angle [deg]: 
        Start=90 # choose fix angle along holder axis
         Stop=0  # do not vary the angle during wobbling
         Step=0  # do not vary the angle during wobbling
    Amplitude [DAC]:
        Start=0.0
         Stop=25.0
         Step=1.0
 press Wobble Button
 ```
 
<details>
  <summary>Wobbler Calibration Output Example</summary>
       <div class="content">
  <img src="https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/wobbler_calib2.png">
</div>
  Standard derivative of wobbler calibration image intensities for number of samples n=[2,3,4] plotted vs maximum tilt magnitude $phi_{max}$. The distance of the maxima scales with n∙2π, therefore the beam tilt amplitude for an arbitrary phase shift can be measured by the distance in the wobbling curve.
    
</details>

---
# Data Evaluation

 ## Specimen drift correction
 
<div class="content">
  <img align="right" src="https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/DriftCorrUI.png">
</div>
The main tool for the specimen drift correction is the digital micrograph script [DriftCorrUI.s]().
Its purpose is to create a custom mask in reciprocal space. This custom mask can constist of combinations of spot- and radial masks to focus the drift correction of the relevant spatial frequencies. After the mask is constructed iwth the user interface it can be applied to the image stack and the drift calculated via cross correaltion or phase correlation function. The drift vector is returned in form of an image, that can be applied to stacks. <br>

The specimen drift correction scheme is proposed in the [publication]() and an [example workflow](#example-data-and-workflow) is given below. 
<br>

<details>
 <summary>Input</summary>
 
 <kbd>Assign</kbd>: Assigns the frontwindow image stack to the Userinterface.<br>
 <kbd>Hanning</kbd>: Apllies a realspace hanning window with radius=0.5* image width/height to the stack center.<br>

</details>

<details>
 <summary>Aperture</summary>
 
 <kbd>Load</kbd>: Assigns frontmost image as working aperture.<br>
 <kbd>Save</kbd>: Saves current working aperture as image.<br>

</details>

<details>
 <summary>Radial masks</summary>
 
 <kbd>Top-head</kbd>: Applies an top-head mask to working aperture. The radius is given in pixel by <kbd>Radius</kbd><br>
 <kbd>B\`worth</kbd>: Applies a Buttherworth to working aperture. The radius is given in pixel by <kbd>Radius</kbd> and its order is defined<br>&emsp;&emsp;&emsp;&emsp; by an user dialog.<br>
 <kbd>Radius</kbd>: . Radius in pixel used for <kbd>Top-head</kbd> and <kbd>B`worth</kbd><br>
 <kbd>Positive</kbd>: Adds a positive spot mask at the most intense pixel inside the current ROI. The radius is given via user<br>&emsp;&emsp;&emsp;&emsp; dialog.<br>
 <kbd>Negative</kbd>: Adds a suppresses FFT signal inside the spot mask at the most intense pixel inside the current ROI.<br>&emsp;&emsp;&emsp;&emsp; The radius is given via user dialog.<br>
 <kbd>Blur</kbd>: Adds Gausssian blur to the edge of an exsiting <kbd>Top-head</kbd> The blur width is given in pixel as user dialog input.<br>
 <kbd>Undo</kbd>: Reverts the last action affecting the working apertue.<br>
    
</details>

<details>
 <summary>Holography</summary>
 
 <kbd>Sideband</kbd>: Removes the most intense bragg peaks inside the sideband. The number of peaks to remove is given as user input dialog.<br>
 <kbd>Fresnel</kbd>: Not fully implemented yet. Should use a gausssian line filter to mask the fresel streak.<br>
    
</details>


<details>
 <summary>Measure</summary>
 
 <kbd>XCF</kbd>: .<br>
 <kbd>PCF</kbd>: .<br><br>
 [&check;] Pairwise: <br><br> 
  <kbd>Max X</kbd>: .<br>
  <kbd>Max Y</kbd>: .<br>
    
</details>

<details>
 <summary>Drift Sequence</summary>
 
 <kbd>Load</kbd>: .<br>
 <kbd>Apply</kbd>: .<br>
    
</details>      
  

 ## Phase shifting series reconstruction

This notation strictly follows the original phase shifting holography reconstruction approach by [Ru et al.](https://doi.org/10.1016/0304-3991(94)90171-6)
Once a phase shifting hologram series of $n$ images has been aqcuired and the [specimen drift correction](#specimen-drift-correction) has been done, the complex exitwave can be reconstructed from the series $I(n)$. $C_i$ describe complex images, that should be recovered from I(n) by solving the least squared problem in matrix form. 
 
 $$ I(n) = C_1 + C_2 \exp \[+i \phi_n\] + C_3 \exp \[ -i \phi_n\],$$
 
 This equation can be rewritten in matrix form:
 
 $$(1,exp[ +i \phi_n],exp[-i \phi_n] ) \left(\begin{array}{c} 
C_1\\
C_2 \\
C_3
\end{array}\right) = I(n) $$ 
 
 Multiplying both sides of the equation by $(1, exp[-i\phi_n], exp[+i \phi_n])^T$, and then
summing both sides of the equation over $n$ leads to a $3×3$ matrix expression:

$$\left( \begin{array}{ccc} N & \sum_{n}\exp(+i \phi_n) & \sum_{n}\exp(-i \phi_n)\\
\sum_{n} \exp(-i \phi_n) & N & \sum_{n}\exp(-2i \phi_n) \\
\sum_{n} \exp(+i \phi_n) & \sum_{n}\exp(+2i \phi_n) & N\\
\end{array} \right)  \left(\begin{array}{c} 
C_1\\
C_2 \\
C_3
\end{array}\right) = \left(\begin{array}{c} 
\sum_{n} I(n)\\
\sum_{n} \exp(-i \phi_n) I(n) \\
\sum_{n} \exp(+i \phi_n)I(n)
\end{array}\right) $$ 

The equation can be solved effeciently by  inverting the $3×3$ matrix. This is done in the [reconstruction Matrix](#reconstruction-matrix) section via a digital micrograph script with an embedded python routine.  

The complex inverted matrix is needed as input to solve the above matrix equation within the [reconstruction Image](#reconstruction-image) step.
 
  ### Measure the carrier frequency phase
  [TiltSeriesUI.s]()<br>
  Input: <li><ul>Phase shifting hologram series as 3D image cube</li></ul>
  Output: <ul><li>1D image containing the $\phi_n$ of each hologram</li> </ul>
  
 The phase shift $\phi_n$ of the hologram carrier frequency of each hologram has to be measured by the <kbd>Measure stack</kbd> of [TiltSeriesUI.s](). The resulting line profile will be used as input for the next step. 
  
  
  ### Reconstruction Matrix
  [ReconPS_Holo_Matrix.s]()<br>
  Input: <ul>
  <li>1D image containing the $\phi_n$ of each hologram</li>
</ul>
  Output: <ul>
  <li>Two 3x3 pixel images containing the real- and imaginary part of the corresponding inverted matrix</li></ul>
  
  ### Reconstruction Image
  [ReconPS_Holo_Image.s]()<br>
  Input: <ul>
<li>1D image containing the $\phi_n$ of each hologram</li><br>
<li>Drift corrected 3D image cube containing the phase shifting series</li><br>
<li>Two 3x3 pixel images containing the real- and imaginary part of the corresponding inverted matrix</li>
</ul>
Output:<ul>

</ul>
  
  ## Reference correction
 
  
  ---
# Example Data and Workflow
---
# Copyright and Disclaimer

<details>
  <summary>Disclaimer</summary>
    
All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.

temscript is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENCE.txt file for any details.

All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.
</details>   
 
