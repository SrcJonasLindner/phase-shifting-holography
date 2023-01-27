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
 <kbd>Stepsize</kbd>:Used stepsize if ramping the angle. Set to zeroi for fixed angle measurements.<br>
  
</details>

<details>
  <summary>Amplitude [DAC]<a name="deg"></a></summary>
    
  <kbd>Start</kbd>: Initial tilt amplitude in units to the digital to analog converter.  <br>
 <kbd>Stop/Its.</kbd>:Final tilt magnitude. For typical biprism voltages from 200-300 V [0,40] is a good first wobbbling<br> &emsp;&emsp;&emsp;&emsp; guess.<br>
 <kbd>Stepsize</kbd>:If stepsize=0: using the stop value as number of measurements to be done. If stepsize>=0 use it to<br>&emsp;&emsp;&emsp;&emsp; determine the steps between <kbd>Start</kbd>and <kbd>Stop</kbd>.<br>
  
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
    
<kbd>Drift</kbd>: Uses the [phase correlation function](https://doi.org/10.1016/S0304-3991(02)00071-2) to roughly calculate the specimen drift  of the frontwindow <br>&emsp;&emsp;&emsp;&nbsp;image stack.<br>
<kbd>Measure stack</kbd>: <br>
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


# Data Evaluation

 ## Specimen drift correction

 ## Phase shifting series reconstruction
  ### Measure the carrier frequency phase
  ### Reconstruction Matrix
  ### Reconstruction Image

  ## Reference correction
  
  ---
# Example Data and Workflow

# Copyright and Disclaimer

All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.

temscript is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENCE.txt file for any details.

All product and company names are trademarks or registered trademarks of their respective holders. Use of them does not imply any affiliation with or endorsement by them.

