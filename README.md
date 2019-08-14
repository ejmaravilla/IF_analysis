# IF_analysis
counting IF cells (Shen Lab)

HALI
#############################
#####::Project Summary::#####

HALI = Hyperoxide Acute Lung Injury

This project will determine the time at which the maximal number of 
alveolar (AT2) cells, putative stem cells, are proliferating after
acute lung injury.

We are inspecting lung tissue sections cut from formalin fixed and 
paraffin embedded blocks of whole lung lobes 
(S=superior, M=medial, I=Inferior, C=cardiac, L=left).

We are marking alveolar cells that are positive for 
surfactant protein C (Sftpc) as putative AT2 cells.

We are marking proliferating alveolar cells positive for Ki67. 

Started: 2019-08-06
Maintainers: Rui Xi, Preetish Kadur L Murthy
Assistant: Erick Maravilla

#########################
#####::Data Origin::#####

This data is acquired at Duke University.
Microscope: Keyence BZ-X series 710
		Located in Dr. You's Lab, CIEMAS 2nd floor
Acquisition Settings:
	Objective: 40x
	Image Area: 2.899 x 2.174 mm (11x11)
	Four channels: BF, Red, Green, Blue
	Per Scan (2.899x2.174mm):
		Time: ~75 minutes
		Images: (11x11) x 4channels = 484 (OR) 121 x 4 = 484
	One image per lobe: - user selected, try to get area with mininal artifacts
	Possibly multiple lobes imaged per slide
	See HALI immunofluorescence protocol for preparation details
