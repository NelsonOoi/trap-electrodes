from trapsim import *
from scipy.integrate import solve_ivp
from scipy.fft import fft, fftfreq

'''
-----------------------------------------------------------------------------------------------------------------------
Initial test from 
https://medium.com/@bldevries/simply-solving-differential-equations-using-python-scipy-and-solve-ivp-f6185da2572d
-----------------------------------------------------------------------------------------------------------------------
'''
# def fun(t, u):
#     x, v = u
#     return [v, -1 * (2*np.pi)**2 *x]

# u_0 = [1, 0]
# N = 10000
# t_max = 15
# t_spacing = t_max/N
# t_span = N * t_spacing

# t_pts = np.linspace(0, t_max, N, endpoint=False)
# relerr = 1.e-20
# abserr = 1.e-20
# result = solve_ivp(fun, (0, 20*np.pi), u_0, t_eval=t_pts,
#                     rtol=relerr, atol=abserr)

# plt.plot(t_pts, result.y[0,:], label = "Numerical solution")
# plt.plot(t_pts, [np.cos(2 * np.pi * t) for t in t_pts], "o", label="Analytical solution")
# plt.xlabel("t")
# plt.ylabel("x(t)")
# plt.xlim(0,5)
# plt.legend()

# yf = fft(result.y[0,:])
# xf = fftfreq(N, t_spacing)[:N//2]
# plt.figure()
# plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), label = "FFT")
# # plt.plot([math.sin(t) for t in t_pts], "o", label="Analytical solution")
# plt.xlabel("f (Hz)")
# plt.ylabel("A")
# plt.xlim(0,10)
# plt.legend()
# plt.grid()

# plt.show()
'''
-----------------------------------------------------------------------------------------------------------------------
New test inspired by Luke's thesis.
-----------------------------------------------------------------------------------------------------------------------
dv/dt = - q/m * d/dx(phi(x))
dx/dt = v
'''

s, electrodes, electrodes_dict = load_trap()
'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 20. # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 35e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)
x_guess = [0., 0., 50.]
x0 = s.minimum(x_guess)
ion_height = x0[2]

files = ['Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv', 'Apr22_gds_Vs_axial.csv', 'Apr22_7dof_gds_Vs_axial.csv']
overall = read_electrode_voltages(files=files)
voltage_set = overall[0]

with s.with_voltages(dcs=voltage_set):
    rf_scale = s.rf_scale(m, q, l, o)
    mu, b = s.mathieu(x0, scale=rf_scale, r=4, sorted=True)
    freqs = mu[:3].imag*o/(2*np.pi)
    modes = b[len(b)//2 - 3:len(b)//2, :3].real
    for i in range(len(freqs)):
        print(f"Mathieu frequency {i+1}:", freqs[i]/1e6, "MHz", modes[i])

    def fun(t, u):
        x, v = u
        dxdt = v
        potential_x = s.potential([x, 0, ion_height], 1).flatten()[0]
        dvdt = - (q / m) * potential_x
        # print('potential', potential_x, 'x', x, 'dvdt', dvdt, 'v', v)
        return [dxdt, dvdt]

    u_0 = [1, 0]
    N = 30000000
    t_spacing = 1e-7
    tmax = N * t_spacing
    t_pts = np.linspace(0, tmax, N)
    t_span = (0, 3*tmax)
    relerr = 1.e-20
    abserr = 1.e-20
    result = solve_ivp(fun, t_span=t_span, y0=u_0, t_eval=t_pts,
                    rtol=relerr, atol=abserr)

    plt.figure()
    plt.plot(t_pts, result.y[0,:] / micron_to_m, label = "Numerical solution")
    # plt.plot([math.sin(t) for t in t_pts], "o", label="Analytical solution")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    # plt.xlim(0,40)
    # plt.legend()
    plt.grid()
    # plt.show()

    yf = fft(result.y[0,:])
    xf = fftfreq(N, t_spacing)[:N//2]
    plt.figure()
    plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), label = "FFT")
    # plt.plot([math.sin(t) for t in t_pts], "o", label="Analytical solution")
    plt.xlabel("f (Hz)")
    plt.ylabel("A")
    plt.xlim(0,5)
    plt.legend()
    plt.grid()
    plt.show()
