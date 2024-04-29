from trapsim import *

'''
Algorithm:

1. Select:
    - interpolation (sin, bezier, param).
    - shuttling time (t_shuttle) between start and finish points.
    - final axial phonon number (n_bar) after cooling runs.
    - axial frequency / well curvature during shuttling. usually keep this constant.
2. Calculate amplitude of initial axial oscillation using nbar * hbar * omega_ax(t) = 1/2 * m_Ca+ * omega_ax(t)^2 * A^2
3. Solve for ion motion over the complete time interval.
   At each time t:
    - Compute position of the well (x_well) using the interpolation functions.
    - Compute axial voltage set for creating a well at x_well.
    - Simulate the axial voltage set to find d/dx(Phi) at the ion's position (x_ion).
4. At the final position, FFT the ion's motion and plot it. Deduce the final axial phonon number from this.
5. ???
6. Profit!

Soon:
- Expand to simulations of multi-ion chains.
- Profile optimizer routine to minimize the final oscillation amplitude.


'''