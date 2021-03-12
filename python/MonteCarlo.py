import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append(".")
from python.Random import Random

# Integrate under the curve sin(pi*x) + 1 from 0 to 1 with all three integration methods:
# Simpson's Rule, Gauss-Legendre quadrature, and Monte Carlo integration
# the analytic solution will be (2 + pi) / pi

if __name__ == "__main__":

	#set default number of steps
	Nsteps = 100

	# set default number of sub-intervals
	Nint = 2

	method = 0

	# read the user-provided arguments from the command line (if there)
	if '-Nsteps' in sys.argv:
		p = sys.argv.index('-Nsteps')
		Nsteps = int(sys.argv[p+1])

	if '-h' in sys.argv or '--help' in sys.argv:
		print ("Usage: %s -Nsteps [number]" % sys.argv[0])
		sys.exit(1)
	
	# just want to automatically compare all three integration methods. only interested in setting the number of steps.

	# function to integrate over, sin(pi*x) + 1
	def fn(x):
		return np.sin(np.pi*x) + 1

	# coefficients for the integration methods
	coeff_NC = [[1.0, 1.0], 
				[1.0, 4.0, 1.0],
				[1.0, 3.0, 3.0, 1.0],
				[7.0, 32.0, 12.0, 32.0, 7.0]]

	# weights, roots taken from Wiki on Gaussian quadrature - https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
	weights_gauss = [[2.0],
					 [1.0, 1.0],
					 [0.555556, 0.888889, 0.555556],
					 [0.347855, 0.652145, 0.652145, 0.347855],
					 [0.236927, 0.478629, 0.568889, 0.478629, 0.236927]]

	roots_gauss = [[0.0],
				   [-0.57735, 0.57735],
				   [-0.774597, 0.0, 0.774597],
				   [-0.861136, -0.339981, 0.339981, 0.861136],
				   [-0.90618, -0.538469, 0.0, 0.538469, 0.90618]]

	# simpson's rule
	def simpson(a, b, n):
		h = (b-a)/n
		prefactor = [h/2.0, h/3.0, 3.0*h/8.0, 2.0*h/45.0][n-1]

		s_int = 0.0
		#print(a,b)
		for i in range(0,n+1):
			xi = a + i*h
			#print("xi, ", xi)
			s_int += prefactor*coeff_NC[n-1][i]*fn(xi)
		return s_int

	# gauss legendre, need to do change of interval - see https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
	def gauss(a, b, n):
		# (b-a / 2) * sum(w_i * f((b-a / 2)*xi + (b+a)/2))
		g_int = 0.0
		for i in range(1, n+1):
			g_int += ((b-a)/2.0) * weights_gauss[n-1][i-1] * fn(((b-a)/2.0)*roots_gauss[n-1][i-1] + ((b+a)/2))
		return g_int

	# monte carlo rejection sampling
	def mc_reject(a, b, x):

		top = np.amax(fn(x)) # get max of function to integrate

		nAccept = 0
		nTotal = 0
		
		# accepted values
		Xaccept = []
		Yaccept = []

		# reject values
		Xreject = []
		Yreject = []

		# sample number
		isample = []
		# calculated values of the integral
		calcInt = []

		random = Random()

		idraw = max(1,int(Nsteps)/100000)
		#print(idraw)
		for i in range(0,Nsteps):
			X = random.rand() # leave as is since x range is [0,1]
			Y = top*random.rand() # transform this to sample from the max value of the function
			#print(top)
			nTotal += 1
			if(Y < np.sin(np.pi*X) + 1.0): #accept if inside
				nAccept += 1
				if(i % idraw == 0):
					Xaccept.append(X)
					Yaccept.append(Y)
			else: # reject if outside
				if(i % idraw == 0):
					Xreject.append(X)
					Yreject.append(Y)
			if(i % idraw == 0):
				isample.append(nTotal)
				calcInt.append(top*nAccept/nTotal)

		return(isample, calcInt)



	# get intervals on [0, 1]
	intervals = []
	rng = np.linspace(0, 1, Nsteps+1)
	for i in range(Nsteps):
		intervals.append((rng[i], rng[i+1]))

	# calculate integrals
	"""integral = 0.0

	if method == 0:
		for interval in intervals:
			integral += simpson(interval[0], interval[1], Nint)
	elif method == 1:
		for interval in intervals:
			integral += gauss(interval[0], interval[1], Nint)
	else:
		print("Something broke!")
		sys.exit(1)"""


	analytic = (2.0 + np.pi)/np.pi

	"""print("Integrating e^x between 0 and 1.")
	print("Analytic answer: ", analytic)
	print("Numerical answer: ", integral)
	print("Analytic-Numerical: ", analytic-integral)"""



	# code to do some systematic comparisons
	
	simp = np.zeros(4)
	gaus = np.zeros(5)

	# get simpson as function of Nint
	for Nint in [1,2,3,4]:
		for interval in intervals:
			simp[Nint-1] += simpson(interval[0], interval[1], Nint)

	# get gauss as function of Nint
	for Nint in [1,2,3,4,5]:
		for interval in intervals:
			gaus[Nint-1] += gauss(interval[0], interval[1], Nint)

	# get monte carlo as function of Nsteps
	mcx = np.linspace(0, 1, Nsteps)
	isample, mc_int = mc_reject(0, 1, mcx)

	print("Comparing Analytic to Monte Carlo Rejection")
	print("Analytic: ", str(analytic))
	print("Monte Carlo: ", str(mc_int[-1]))
	print("Simpson (n=2): ", str(simp[1]))
	print("Gauss-Legendre (n=2): ", str(gaus[1]))
	print("Analytic - Monte Carlo: ", str(analytic - mc_int[-1]))
	print("Analytic - Simpson: ", str(analytic - simp[1]))
	print("Analytic - Gauss-Legendre: ", str(analytic - gaus[1]))

	#plot calculated pi vs sample number
	fig1 = plt.figure()
	plt.plot(isample, mc_int)
	plt.ylabel(r'Approximate integral')
	plt.xlabel("Sample number")
	plt.xlim(0,isample[len(isample)-1])
	ax = plt.gca()
	ax.axhline(y=analytic,color='green',label=r'true integral $\sim {:.6f}$'.format(analytic))
	ax.axhline(y=simp[2],color='red',label=r'Simpson (n=2) $\sim {:.6f}$'.format(simp[2]))
	ax.axhline(y=gaus[2],color='blue',label=r'Gauss (n=2) $\sim {:.6f}$'.format(gaus[2]))
	plt.title(r'Approximation of $\int_0^1 \! sin(\pi x) + 1 \mathrm{d}x$ as a function of number of samples')
	plt.legend()
	plt.show()
	fig1.savefig("compare_ints.jpg", dpi=180)


# DISCUSSION:
# For a relatively simple function as sin(pi*x) + 1, the deterministic methods perform much better at even modest order (n=2) than the monte-carlo method.
# The monte-carlo method converges on the correct answer, but with much higher error than the deterministic methods, and I don't think it will reach that precision
# in anything remotely resembling reasonable time (for this specific integral)