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
<div class="content">
  <img align="right" src="https://github.com/SrcJonasLindner/phase-shifting-holography/blob/main/doc_images/TiltSeriesUI.png">
</div>
  
  
  
  ## Calibrating the tilt
    
    The direction perpendicular to the biprism orientation is the most sensitive tilt direction. To minimize the coupling from mechanical vibrations into the TEM holder, it is (often) favoarble for the coherence to choose the biprism orientation perpendictlar to the long holder axis. In this special case,the x-direction of the gun tilt lens points along the biprsim. In <mark>TiltSeriesUI.s<mark> the an angle of 0 degree corresponds to the holder axis, so 90 degree is the favorable tilt direction, but the UI allows to tilt in any particular direction, depending on the biprism settings used.
    
    If one tilt angle has been chosen, the magnitude of the tilt has to be calibrated to the used biprism voltage. In our paper we use a medium voltage of 250V, which was sufficient  to seperate the frist order bragg reflections of center- and sideband of platinum, while providing a good visibility of 20%. Goal of the calibration process is to find the beam tilt magnitude in the particular direction, that corresponds to a phase shift of 2PI in the hologram. 
    
    
    
# Data Evaluation

 ## Specimen drift correction

 ## Phase shifting series reconstruction
  ### Measure the carrier frequency phase
  ### Reconstruction Matrix
  ### Reconstruction Image

  ## Reference correction
  
  ---
# Example Data and Workflow
