# How to Download
To use the MATLAB GUI, you'll need to install *parametricpropulsionmodel/PTPM.mlappinstall.* Opening this file while in MATLAB will install the GUI as an app and will be ready to use.

# Background
This project aims to develop a tool for modeling turbofan performance for aircraft conceptual designers. This will benefit aircraft designers as current literature and available resources on propulsion performance are limited due to proprietary restrictions. The model developed is specific to subsonic two-spool turbofans with separate exhausts and no afterburners. Objectives of our research are as follows:
- Utilize easily accessible published manufacturer data
- Return approximate thrust and specific fuel consumption trends
- Create an easy to use tool for any designer knowledge
- Allow future expansion of model

The two important desired outputs of the model are thrust and specific fuel consumption (SFC). For a particular engine, knowing how the thrust and SFC vary with mach number and altitude is a necessity for aircraft conceptual designers. Even low-fidelity models will enable engineers to choose an initial engine for an aircraft and model the performance of the aircraft on a specified mission.  

# Methodology
Software development began with research on existing models. Models exist for a variety of engine types, many with varied assumptions and accuracies, but exist mostly as a preliminary tool for an engine designer, not an aircraft designer. Models, typical engine parameter ranges, relationships between parameters, and underlying assumptions were studied until model development began.​

Model development begins with the on-design condition. Designers will design to the most critical condition and display associated component parameters in engine brochures. On-design analysis is developed at sea-level static conditions to obtain information about engine component performance. The primary purpose of developing an on-design model is to solve for component efficiencies to use in the following off-design calculations.  ​

Applying component information and several key assumptions to component adjustment to flight conditions, the off-design analysis is developed to model how each component changes under a different environment.  With this analysis, engine performance can be modeled across the operational envelope. ​

Air flow through the engine is modeled with specific heat varying with temperature and fuel to air ratio, which in turn influences thermodynamic properties throughout the engine. A driving algorithm used throughout our analysis, known as FAIR, models this variation with NASA thermochemical data and is used to find all thermodynamic properties based simply on the fuel to air ratio and either 1) temperature, 2) pressure, or 3) enthalpy.​

By implementing an on-design analysis once and an off-design analysis at a specified range of mach numbers and altitudes, the desired thrust and SFC data for the range of operating conditions can be obtained.

# Discussion and Conclusions
Minimal access to engine data was a particular hindrance throughout the model development. Several partial additional data sets were found, but none that could grant full confidence in the robustness and accuracy of the software. Future work should include testing accuracy with additional engine data and map out error in the analysis. The results from the DGEN-390 indicate accuracy is acceptable for a conceptual design, however additional confirmation should occur to prove reliable for actual aircraft designers.

In order for the method to be easily used, a graphical user interface (GUI) was designed. The GUI allows users of all knowledge levels to obtain meaningful data from just 9 inputs of easy to obtain manufacturer published data. Each input is annotated with helpful information such as parameter definitions, realistic values, and effects on the system and analysis. Parameters are primarily obtained from manufacturer brochure data, but some estimations must be made by the user. If the user is unsure of a value, the user can leave it blank and the software will guess and use a realistic value.
​
Maximum thrust and specific fuel consumption data is displayed in a visual plot. The curves are fully customizable over a range of Mach numbers and altitudes. Results are also tabulated. More detailed engine data such as thermodynamic states across engine station and component performance parameters is stored by the program and can be accessed via the MATLAB workspace.​ Data can also be exported in custom .mat or .txt files so engine data can be used in other aircraft design analysis methods.
​
It is also evident that some combinations of inputs result in a nonfunctional engine. For this reason, multiple error codes are in place to allow the user to troubleshoot invalid inputs. Fidelity of the analysis can be adjusted to vary the number of data points, and the subsequent runtime to user needs.

The development of this engine analysis software will benefit the aircraft conceptual design process by providing initial low-fidelity trends of thrust and specific fuel consumption to engineers. One of the most groundbreaking aspects of the software is eliminating the need for proprietary engine information to obtain accurate trends. In addition, the software is easy to use and provides accurate data with minimal inputs. For these reasons, the parametric turbofan propulsion model will benefit aircraft conceptual designers greatly.


# More Information
https://surp.calpoly.edu/2020/parametric-turbofan-propulsion-model/
