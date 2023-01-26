# Phase-shifting-holography

Plaese read the Ultramicroscopy publication < paper doi>

Software package collecting scripts for the data collection and evaluation of phase shifting holography image stacks with Digital Micrograph Scripting and TitanScripting for electron microscopy.

# Requirements
  ## Digital Micrograph 3.5X+ with Python support installed
  
  Install the latest version DM-Version from: [https://www.gatan.com/installation-instructions]
  
  ## Software Microscope control to lenses (optional, but strongly recommended)
  The phase shifting measurement procedure involes a routine to use the GunTilt Lens of the Titan to tilt the beam off from the optical axis to create phase    shifts in the holograms, that are the heart of this method. There is a default option using the Digital Mricorgraph Microscope interface, but with our Microscope this procedure was way too slow and needed approx 1-2 seconds to create one of the beam tilts. Since a phase shifting holo series covers around several tens of images, the specimen driftin the dead time is a serious problem. 
  
  
  Tem scripting interface (https://temscript.readthedocs.io/) or the ability to precompile microscope commands as executable files, that then are included in the measurement script.
  

# Data collection

# Data Evaluation

# Example Workflow
