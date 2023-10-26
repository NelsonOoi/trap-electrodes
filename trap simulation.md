Oct 25


Discussion Oct 23, 2023
- Electrode potential fitting quadratic parameters
- Vighnesh suggests creating two perpendicular planes to do the fit
	- Coordinate axis: z trap axis, y above trap plane, x perpendicular to trap axis
	- Use xy plane and xz plane to do fitting
	- In each plane the third parameter is set to a constant


Progress log 1
Oct 20
- Working on applying unit voltages on each electrode and simulating potential
- Important files
	- system.py - defines System class that one can add electrodes to
		- individual_potential seems to do what we need it to do, by applying a unit voltage to each
		- electrical_potential
- from OneNote: In the direction out of the page (z axis), the ion is at z=56.5 µm above the waveguide plane (this is 50 µm above the top metal, and the chip stack is approximately 6.5 µm thick above the waveguides).

Sept 30
- Created function to generate trap electrodes
- Referred to 2023_09_18_Quetzal_Reticle_Processed_wLabels.gds file for electrode dimensions
	- DC electrode widths: 120 µm
	- DC electrode height: 1000 µm
	- Mid DC electrode height: 20 µm
	- Mid RF electrode height: 70 µm
	- Spacings: 5 µm
	- Effective widths are spacings/2 + widths

![[quetzal.png]]