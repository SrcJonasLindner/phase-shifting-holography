# Phase-shifting-holography

Plaese read the Ultramicroscopy publication < paper doi>

Software package collecting scripts for the data collection and evaluation of phase shifting holography image stacks with Digital Micrograph Scripting and TitanScripting for electron microscopy.

--- 
** Table of contents: **

1. [Requiements](#requirements)
1. [Data Collection](#data-collection)
  1.1 Tilt Series UI(#Tilt-series-ui)
  1.1 Calibrating the tilt(#Calibrating-the-tilt)
1. [Data Evaluation](#data-Evaluation)
 1.1 [Specimen drift correction](#specimen-drift-correction)
 1.1 [Phase shifting series reconstruction](#phase-shifting-series-reconstruction)
    1.1.1 [Measure the carrier frequency phase](#measure-the-carrier-frequency-phase)
    1.1.1 [Reoncstruction Matrix](#reconstruction-matrix)
    1.1.1 [Reconstruction Image](#reconstruction-image)
1. [Example Data and Workflow](#Example-Data-and-Workflow)

# Requirements
  ## Digital Micrograph 3.5X+ with Python support installed
  
  Install the latest version DM-Version from: [https://www.gatan.com/installation-instructions]
  
  ## Software Microscope control to lenses (optional, but strongly recommended)
  The phase shifting measurement procedure involes a routine to use the GunTilt Lens of the Titan to tilt the beam off from the optical axis to create phase    shifts in the holograms, that are the heart of this method. There is a default option using the Digital Mricorgraph Microscope interface, but with our Microscope this procedure was way too slow and needed approx 1-2 seconds to create one of the beam tilts. Since a phase shifting holo series covers around several tens of images, the specimen driftin the dead time is a serious problem. 
  
  
  Tem scripting interface (https://temscript.readthedocs.io/) or the ability to precompile microscope commands as executable files, that then are included in the measurement script.
  

# Data collection
  ## Tilt Series UI
  ## Calibrating the tilt
# Data Evaluation

 ## Specimen drift correction

 ## Phase shifting series reconstruction
  ### Measure the carrier frequency phase
  ### Reoncstruction Matrix
  ### Reconstruction Image

  ## Reference correction
  
  
# Example Data and Workflow
