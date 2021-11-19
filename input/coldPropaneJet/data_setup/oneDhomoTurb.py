#####################################################################
# DOL 10/17/11 (original in Matlab; converted to Python 2/24/16)
#
# Initialize a 1-D turbulent velocity field
# Based on Pope's model spectrum.  See "Turbulent Flows" p. 232
# Compute the parameters cl_p and ce_p in Pope's model turbulent
#    kinetic energy spectrum.  See Pope "Turbulent Flows" p. 232
# These parameters are enforced by integrating the spectrum to
#    get the desired energy and dissipation rate.
#
# IFFT: u_j    =       sum n=0,N-1 uhat_n * exp( 2*pi*i*n*j/N)
# FFT:  uhat_n = 1/N * sum j=0,N-1 u_j    * exp(-2*pi*i*n*j/N)
# If we have N points, then we get N/2+1 waves for even N.  See Fig. 2 of the dft.pdf writeup.
# Parseval's theorm says that [1/N * sum uhat*conj(uhat)] = [sum u*u], sums are 0,N-1
# We have: (sumh(E*dk)) = (1/2*u'u') = (1/2/N*sum(u*u)) = (1/2/N/N*sum(uhat*conj(uhat))) =
#          (1/N/N*sumh(uhat*conj(uhat))), where sumh means sum from 0 to N/2 (h for half).
# Then the first and last equalities are (sumh(E*dk)) = (1/N/N*sumh(uhat*conj(uhat))).
# To get u, we need uhat so we can take u=ifft(uhat).  uhat comes from E, and they are
#   related by the above equality.  We then take the amplitude of uhat as A = sqrt(E*dk*N*N).
# We have N/2+1 of these.
# We then randomize the phase angle phi.
# We then compute uhat = a + bi using a=A*cos(phi) and b=A*sin(phi).
# Then we use conjugate symmetry to get uhat for n=N/2+1 to N-1.
# Then we get u.
#
#####################################################################

from __future__ import division
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad, quadrature, romberg
import matplotlib.pyplot as plt

class homoturb :

    def __init__(self, nu=1.5E-5, rho=1.2, L=1.0, Li=0.1, up=6.65, eta=0.001, specifyUprime=False ) :

        s = self

        #----------------------- User specifications

        s.nu   = nu                    # kinematic viscosity
        s.rho  = rho                   # density
        s.L    = L                     # domain size
        s.Li   = Li                    # integral scale (specifies ke)
        s.neta = 2                     # number of points per eta scale (integer)

        #------ set one of these and the flag specifyUprime

        s.up  = up                     # velocity fluctuation (specifies ke)
        s.eta = eta                    # Kolmogorov length scale
        ### specifyUprime: true uses uprime, false uses eta

        #----------------------- Useful values

        if specifyUprime :
            s.kine = 1/2*s.up**2
            s.eps  = s.kine**(3/2)/s.Li
            s.eta  = (s.nu**3/s.eps)**0.25      # Kolmogorov length
        else :
            s.eps  = s.nu**3/s.eta**4
            s.kine = (s.Li*s.eps)**(2/3)
            s.up   = np.sqrt(2*s.kine)

        s.Re      = (s.Li/s.eta)**(4/3)
        s.ueta    = (s.eps*s.nu)**0.25
        s.lamb    = np.sqrt(10)*s.eta**(2/3)*s.Li**(1/3)
        s.Re_lamb = s.up*s.lamb/s.nu
        s.mu      = s.nu*s.rho
        s.tau_i   = s.Li/s.up
        s.tau_eta = s.eta/s.ueta

    #----------------------- set the energy spectrum

    def setEspectrum(self) :

        s = self

        s.setSpectrumParameters()

        #------------------------ Set the grid

        s.N  = int(np.ceil(s.neta*s.L/s.eta))
        s.dx = s.L/s.N                        # physical grid spacing
        s.x  = np.arange(0,s.N)*s.dx          # physical grid

        s.k  = np.arange(0,s.N)/s.L           # wave grid
        s.dk = s.k[1]-s.k[0]                  # wave grid spacing
        s.n  = np.arange(0,s.N)               # wave space indicies
        s.nh = np.arange(0,int(s.N/2+1))      # half of the wave space indicies
        s.Nh = int(np.floor(s.N/2+1))         # half the grid points: N=6 --> 4, N=5 --> 3
        s.k[0] = 1.0E-50*s.k[1]               # avoid division by zero below

        #----------------------- Set turbulence parameters and spectrum

        fl  = ( s.k[s.nh]*s.Li / np.sqrt((s.k[s.nh]*s.Li)**2 + s.cl_p) )**(5/3+2)
        fe  = np.exp( -5.2 * ( ( (s.k[s.nh]*s.eta)**4 + s.ce_p**4 )**0.25 - s.ce_p ) )
        s.E = 1.5*s.eps**(2/3)*s.k[s.nh]**(-5/3)*fl*fe
        s.E[0] = 0

    #---------------------- Get u profile

    def getUprofile(self) :

        s = self

        #----------------------- Get the fourier coefficients (uhat)

        A   = np.sqrt(s.E*s.dk*s.N*s.N)         # mode amplitude
        phi = np.random.rand(s.Nh)*2*np.pi      # randomize the phase angle

        a = A*np.cos(phi)                       # get the real part
        b = A*np.sin(phi)                       # get the imag part

        s.uhat = np.zeros(s.N,complex)
        s.uhat[s.nh] = a + 1j*b                  # form the complex modes

        #----------------------- Impose conjugate symmetry

        ii = np.arange(1, int(np.ceil(s.N/2)))  # imposing conjugate symmetry
        jj = s.N-ii                             # for N=6: 1234 --> 123432; for N=5: 12332
        s.uhat[jj] = np.conj(s.uhat[ii])

        #------------------------ Grab the u profile

        s.u = np.fft.ifft(s.uhat)                  # get the velocity field
        s.u = np.real(s.u)

    #----------------------- Output some sanity checks and plot

    def output(self, LdoPlot=False) :

        s = self

        print("L          = %g" %(s.L))
        print("Li         = %g" %(s.Li))
        print("up         = %g" %(s.up))
        print("kine       = %g" %(s.kine))
        print("eps        = %g" %(s.eps))
        print("eta        = %g" %(s.eta))
        print("Re         = %g" %(s.Re))
        print("ueta       = %g" %(s.ueta))
        print("lamb       = %g" %(s.lamb))
        print("Re_lamb    = %g" %(s.Re_lamb))
        print("mu         = %g" %(s.mu))
        print("tau_i      = %g" %(s.tau_i  ))
        print("tau_eta    = %g" %(s.tau_eta))
        print("cl_p       = %g" %(s.cl_p))
        print("ce_p       = %g" %(s.ce_p))

        print('\nThe next four values should be the same:')

        print("upup_div_2                = %g" %(s.up*s.up/2))                       # average kinetic energy specified
        print("sum_Edk                   = %g" %(sum(s.E*s.dk)))                     # integral of E = 1/2*up*up
        print("sum_uu_div_2N             = %g" %(sum(s.u*s.u)/2/s.N))                # average kinetic energy result
        print("sum_uhat_conjuhat_div_2NN = %g" %(sum(np.real(s.uhat*np.conj(s.uhat)))/2/s.N/s.N)) # fourier energy
        print("\n\n")

        if LdoPlot :

            plt.subplot(3,1,1)
            plt.loglog(s.k[s.nh]*s.eta, s.E/s.eta/s.ueta**2)
            plt.xlim([0.001,10])
            plt.ylim([0.01,10000])
            plt.xlabel(r'$\kappa\eta$');
            plt.ylabel(r'$E/\eta u_{\eta}^2$');

            plt.subplot(3,1,2)
            plt.plot(s.k,np.abs(s.uhat))
            plt.xlabel('Wavenumber (1/m)');
            plt.ylabel('|uhat|');

            plt.subplot(3,1,3)
            plt.plot(s.x,s.u)
            plt.ylim([np.min(s.u)+np.min(s.u)*0.1, np.max(s.u)*(1.1)])
            plt.xlabel('Position (m)');
            plt.ylabel('Velocity (m/s)');

            plt.tight_layout()
            plt.show()

    #----------------------- Get spectrum parameters

    def setSpectrumParameters(self) :

        s = self

        def fGetParams(cparams) :

            def espec(k) :
                fl = ( k*s.Li / np.sqrt((k*s.Li)**2 + s.cl_p) )**(5/3+2)
                fe = np.exp( -5.2 * ( ( (k*s.eta)**4 + s.ce_p**4 )**0.25 - s.ce_p ) )
                return 1.5*s.eps**(2/3)*k**(-5/3)*fl*fe

            def twonuk2E(k) :
                return 2*s.nu*k*k*espec(k)

            s.cl_p   = cparams[0]
            s.ce_p   = cparams[1]

            int1 = quad(espec,    0, np.inf)[0]
            int2 = quad(twonuk2E, 0, np.inf)[0]

            return np.array([int1-s.kine, int2-s.eps])

        s.cl_p = 3.066   #5.9787              # initial guesses for params
        s.ce_p = 0.457   #0.4036
        cparams = np.array([s.cl_p, s.ce_p])
        cparams = fsolve(fGetParams, cparams)
        s.cl_p, s.ce_p = cparams[0], cparams[1]

    #----------------------- Set turbulence parameters and spectrum


#####################################################################

# h = homoturb()
# h.setEspectrum()
# h.getUprofile()
# h.output(LdoPlot=True)

