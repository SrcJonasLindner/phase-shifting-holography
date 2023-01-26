# Phase-shifting-holography

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

---

# Requirements

  ## Digital Micrograph 3.5X+ with Python support installed
  
  Install the latest version [Digital Micrograph](https://www.gatan.com/installation-instructions) version
  
  ## Software Microscope control to lenses (optional, but strongly recommended)
  The phase shifting measurement procedure involes a routine to use the GunTilt lens of the Titan to tilt the beam off from the optical axis to create phase    shifts in the holograms, that are the heart of this method. There is a [default option](http://www.dmscripting.com/tem_control.html) using the Digital Micrograph Microscope interface, but with our Microscope this procedure was way too slow and needed approx 1-2 seconds to create one of the beam tilts. Since a phase shifting holo series covers around several tens of images, the specimen driftin the dead time is a serious problem. 
 
  [Tem scripting interface](https://temscript.readthedocs.io/) or the ability to precompile microscope commands as executable files, that then are included in the measurement script.
  
---
# Data collection

The acquisition of the phase shifting holo image series requires a calibration of the magnitude and direction of the beam tilt. Both, calibration and the image series recording are implemented in the <mark>TiltSeriesUI.s<mark>. 
---
  ## Tilt Series UI
    
<details>
  <summary>Acquire</summary>
</details>

<details>
  <summary>Angle [deg]</summary>
</details>
    
<details>
  <summary>Amplitude [DAC]</summary>
</details>
    
<details>
  <summary>Options</summary>
</details>
    
<details>
  <summary>Evaluate</summary>
</details>    
    
<div class="content">
  <img align="right" src="https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/TiltSeriesUI.png">
</div>
  
  
  
  ## Calibrating the tilt
    
FFor a given biprism voltage and a certain beam tilt angle the magnitude $\phi$ of the beam tilt in DAC-units has to be calibrated. To minimize the influence of mechanical vibration the biprism should be oriented perpendicular to the holder axis. In the Titan the coordinate system of the beam deflection the y-direction is parallel to the holder axis. The wobbling calibration procedure takes a series of images. The exposure time of each image is synchronized with a tilt-ramp during this exposure. The number of points within the ramp $n$ (samples) has to be chosen beforehand. Each image is acquired with a different maximum tilt $\phi_{max}$. The intensity of the acquired wobbling-image can be described as:
    
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
   
# Data Evaluation

 ## Specimen drift correction

 ## Phase shifting series reconstruction
  ### Measure the carrier frequency phase
  ### Reconstruction Matrix
  ### Reconstruction Image

  ## Reference correction
  
  ---
# Example Data and Workflow
