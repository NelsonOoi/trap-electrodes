Electrode simulation goals
- [x] trap chip gds parser
	- [x] maps electrode names to polygons
- [x] finds curvatures directly from individual electrode potentials
- [x] solves for basis voltage sets
	- axial
	- x,y-comp
	- tilt at user specified angle
- [x] target frequencies
- [x] radial mode breaking degeneracy
- transport waveforms
- export
	- coefficient fits
	- basis voltage sets
	- overall tilt+confining voltage sets

Based on work - Feb 18 [[todo]]
1. Set RF = 30MHz, and find the corresponding $V_{peak}$ that gives $f_{pp,y} \approx f_{pp,z} \approx 3MHz$. (Radial frequencies)
	- Using the freqs_pp calculation, we find that the required curvature for 3MHz is:
		- $u_r = \frac{\omega_r^2 m}{e} \approx \frac{(2\pi \times 3\times10^6)^2(40 * 1.67 \times 10^{-27})}{1.6\times10^{-19}}$
2. Work backward to find the $V_{dc}$ that gives $f_{dc, x} \approx 1MHz$.
	1. Have curvatures corresponding to (+2, -1, -1).
	2. $V_{dc}$ simply scales this curvature.
		1. Need analytical expression for this?
		2. $\omega_x = \sqrt{\frac{eu_x}{m}} \implies u_x = \frac{\omega_x^2 m}{e} \approx \frac{(2\pi \times 10^6)^2(40 * 1.67 \times 10^{-27})}{1.6\times10^{-19}}$
3. Compare 5 & 7 electrodes-a-side.
	1. Does one require lower voltages than the other? Constraint is 10V peak.


Jan 1, 2024
- Currently debugging why different configurations don't result in the expected metastable potential along z-axis.
- Steps so far:
	- Extending measurement axis along z-direction reveals that the turnaround point occurs around 100µm to 300µm for an applied unit voltage.



I utilized the NIST Electrode Python package to generate voltage sets for ion confinement in a surface-electrode ion trap. The results demonstrated that it was possible to obtain reasonable trap depth through techniques described in Allcock et. al.. I also worked with Andrea, my partner on this project, to compare different methods (Python vs COMSOL) for generating analytical solutions for the trap. Concurrently, I set up the beam path for a fiber noise cancelation system using optics equipment and RF signal modulation equipment; I was able to characterize the output signal from the interfered beams.

Discussion Nov 13, 2023
- Fixed bugs in fitting code
	- Was previously using y-axis values to fit z, that's fixed
	- Fixed ```A_grouped``` usage in the linear solver, it should be the transpose instead
- Realizations about Allcock paper
	- 6 voltage groups are used because there are 6 degrees of freedom in the coefficients (x1, x2, y1, y2, z1, z2)
	- This allows a square matrix to be formed and can be solved using a linear solver
- Observations
	- Symmetry in the electrode contributions to potential curvature is present
	- Using these symmetries, the coefficients may cancel out depending on grouping
		- Example: first order x-dependence vanishes with symmetric grouping for axial voltage set 
		- Groupings that result in complete cancelation reduce the degrees of freedom
	- Therefore the ability to use the linear solver depends heavily on which grouping is used
- Ideas
	- Utilize only the 2nd order coefficients for the linear solver when dealing with endcap voltages
		- Does this sacrifice anything by way of neglecting the first order coefficients?
- Questions
	- Why does the z-axis potential slope up and then decrease when moving along the +z direction?
		- Seems that this makes sense as the field lines point into the trap surface when near the surface, and point away as you move up above the surface
	- 

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