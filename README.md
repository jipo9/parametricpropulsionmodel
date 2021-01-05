# 1D Parametric Turbofan Propulsion Model
Our project develops an approximation tool for modeling propulsion for aircraft preliminary design. 
Our model is specific to subsonic two-spool turbofans with separate exhausts and no afterburners. Objectives of our research are as follows:
- Utilize easily accessible published manufacturer data
- Return approximate thrust and specific fuel consumption trends
- Create an easy to use tool for any designer knowledge
- Allow future expansion of model

# How to Download
To use the MATLAB GUI, you'll need to install *parametricpropulsionmodel/GUI/PTPM.mlappinstall.* This should install an app in MATLAB.

# Methodology
Software development began with research on existing models. Models exist for a variety of engine types, many with varied assumptions and accuracies, but exist mostly as a preliminary tool for an engine designer, not an aircraft designer. Models, typical engine parameter ranges, relationships between parameters, and underlying assumptions were studied until model development began.
Model development begins with the on-design condition. Designers will design to the most critical condition and display associated component parameters in engine brochures. On-design analysis is developed at sea-level static conditions to obtain information about engine component performance. The primary purpose of developing an on-design model is to solve for component efficiencies to use in the following off-design calculations.
Applying component information and several key assumptions to component adjustment to flight conditions, the off-design analysis is developed to model how each component changes under a different environment.  With this analysis, engine performance can be modeled across the operational envelop. ​
Air flow through the engine is modeled with specific heat varying with temperature and fuel to air ratio, which in turn influences thermodynamic properties throughout the engine. An algorithm known as FAIR models this variation with NASA thermochemical data and is used throughout the analysis to find all thermodynamic  properties based simply on the fuel to air ratio and either 1)temperature, 2)pressure, or 3)enthalpy.
Finally, the GUI is created to guide a user through the process of using the software. Using MATLAB's App Designer, the app will provide a user-friendly interface to run our analysis. Effort is made to provide the user with appropriate information to allow for informed decisions without overwhelming with extraneous information.​
