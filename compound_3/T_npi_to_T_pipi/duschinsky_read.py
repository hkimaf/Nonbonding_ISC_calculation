import numpy as np
import copy


duschinsky = np.loadtxt('duschinsky_formatted.txt')
displacement = np.loadtxt('displacement_formatted.txt')
freq_initial = np.loadtxt('t_npi_formatted.txt')
freq_final = np.loadtxt('t_pipi_formatted.txt')
d_fi = np.loadtxt('d_fi.txt')
adiabatic = - 0.005789960000016
atomic_time = 2.418884326505e-17
c_cm = 2.998e10
boltzmann = 3.166811563455607e-06
beta = 1 / boltzmann / 300.

partition = 1.0

nmodes = len(displacement)

omega_i = np.array([[0. for i in range(nmodes)] for j in range(nmodes)])
omega_f = copy.copy(omega_i)

R_matrix = np.array([[0. for i in range(nmodes)] for j in range(nmodes)])
for i in range(nmodes):
	for j in range(nmodes):
		R_matrix[i][j] = d_fi[i] * d_fi[j]

for i in range(nmodes):
	partition *= np.exp(- 0.5 * freq_initial[i] * beta) / (1 - np.exp(- freq_initial[i] * beta))
print('partition is', partition)

class IC:
	def __init__(self, nmodes, duschinsky, displacement, freq_init, freq_final, time, beta, previous_sqrt, partition, adiabatic, NAC_matrix):
		self.beta = beta
		self.time = time
		self.nmodes = nmodes
		self.ifreq = freq_initial
		self.ffreq = freq_final
		self.duschinsky = duschinsky
		self.displacement = displacement
		self.previous_sqrt = previous_sqrt
		self.partition = partition
		self.adiabatic = adiabatic
		self.NAC_matrix = NAC_matrix
	def amatrix(self, nmodes, time):
		aplus = np.array([[0. + 0. * 1j for i in range(nmodes)] for j in range(nmodes)])
		aminus = copy.copy(aplus)
		A = copy.copy(aminus)
		for i in range(nmodes):
			tilde_n = 1 / (np.exp(beta * freq_initial[i]) - 1)
			aminus[i][i] = (tilde_n + 1) - tilde_n * np.exp(freq_initial[i] * time * 1j)
			aplus[i][i] = (tilde_n + 1) + tilde_n * np.exp(1j * freq_initial[i] * time)
			A[i][i] = aminus[i][i] / aplus[i][i]
		return (A, aplus, aminus)

	def cmatrix(self, nmodes, time):
		cplus = np.array([[0. + 0. * 1j for i in range(nmodes)] for j in range(nmodes)])
		cminus = copy.copy(cplus)
		for i in range(nmodes):
			cminus[i][i] = 1. - np.exp(- 1j * freq_final[i] * time)
			cplus[i][i] = 1. + np.exp(- 1j * freq_final[i] * time)
		return (cplus, cminus)

	def XI(self, aplus, aminus, cplus, cminus):
		xi = np.dot(np.dot(aplus + aminus, aplus - aminus), cplus - cminus)
		det = np.linalg.det(xi)
		return det
	def L(self, aplus, aminus, omega_f_inv, omega_i, duschinsky_T, duschinsky_inv, cplus, cminus):
		L = np.dot(np.dot(np.dot(np.dot(cminus, omega_f_inv), duschinsky_T), omega_i), aplus) + np.dot(np.dot(cplus, duschinsky_inv), aminus)
		LP = np.dot(np.dot(np.dot(np.dot(cplus, omega_f_inv), duschinsky_T), omega_i), aminus) + np.dot(np.dot(cminus, duschinsky_inv), aplus)
		return (L, LP)
	def J(self, aplus, aminus, cplus, cminus, omega_i_inv, omega_f, duschinsky, duschinsky_Tinv):
		J = np.dot(np.dot(np.dot(np.dot(aplus, omega_i_inv), duschinsky_Tinv), omega_f), cminus) + np.dot(np.dot(aminus, duschinsky), cplus)
		return J
	def LAMBDA(self, aplus, aminus, cplus, cminus, L_inv, LP_inv):
		Lambda = np.dot(np.dot(aplus, L_inv), cplus) - np.dot(np.dot(aminus, LP_inv), cminus)
		return Lambda
	def SQRT(self, det, partition, SS_det, JL_det, sqrt_prev):
		sqrt = np.sqrt(det / partition / (SS_det * JL_det) / partition)
		sqrt_difference = sqrt_prev - sqrt
		sqrt_sum = sqrt_prev + sqrt
		sqrt_difference_abs = np.abs(sqrt_difference)
		sqrt_sum_abs = np.abs(sqrt_sum)
		
		if (sqrt_prev == 0. + 0. * 1j):
			sqrt_updated = sqrt
		else:
			if (sqrt_difference_abs > sqrt_sum_abs):
				sqrt_updated = - sqrt
			else:
				sqrt_updated = sqrt
		return sqrt_updated
	def rho(self):
		A_set = self.amatrix(self.nmodes, self.time)
		A = A_set[0]
		aplus = A_set[1]
		aminus = A_set[2]
		C_set = self.cmatrix(self.nmodes, self.time)
		cplus = C_set[0]
		cminus = C_set[1]
		omega_f = np.array([[0. for i in range(self.nmodes)] for j in range(self.nmodes)])
		omega_i = copy.copy(omega_f)
		omega_f_inv = copy.copy(omega_f)
		omega_i_inv = copy.copy(omega_f)
		duschinsky_T = np.transpose(self.duschinsky)
		duschinsky_inv = np.linalg.inv(self.duschinsky)
		duschinsky_TInv = np.linalg.inv(duschinsky_T)
		for i in range(self.nmodes):
			omega_i[i][i] = self.ifreq[i]
			omega_f[i][i] = self.ffreq[i]
			omega_i_inv[i][i] = 1 / self.ifreq[i]
			omega_f_inv[i][i] = 1 / self.ffreq[i]
		det = self.XI(aplus, aminus, cplus, cminus)
		L_set = self.L(aplus, aminus, omega_f_inv, omega_i, duschinsky_T, duschinsky_inv, cplus, cminus)
		L = L_set[0]
		LP = L_set[1]
		L_inv = np.linalg.inv(L)
		LP_inv = np.linalg.inv(LP)
		J = self.J(aplus, aminus, cplus, cminus, omega_i_inv, omega_f, self.duschinsky, duschinsky_TInv)
		Lambda = self.LAMBDA(aplus, aminus, cplus, cminus, L_inv, LP_inv)
		zeta = np.dot(np.dot(A, omega_i), self.displacement)
		eta = np.dot(np.dot(np.linalg.inv(J), aminus), self.displacement)

		psi = self.displacement - np.dot(np.dot(self.duschinsky, cplus), eta)
		SS = np.dot(duschinsky_T, self.duschinsky)
		JL = np.dot(J, L)
		SS_det = np.linalg.det(SS)
		JL_det = np.linalg.det(JL)
		SIL = np.dot(np.dot(duschinsky_T, omega_i), Lambda)
		FCE = np.dot(np.dot(omega_f, cminus), eta)
		SIAP = np.dot(np.dot(np.dot(duschinsky_T, omega_i), A), psi)
		zeta_dot_psi = np.dot(zeta, psi)
		sqrt_updated = self.SQRT(det, self.partition, SS_det, JL_det, self.previous_sqrt)
		common_factor = (- 1 / np.sqrt(2.)) ** self.nmodes * sqrt_updated * np.exp(- zeta_dot_psi)
		rho = np.array([[0. + 0. * 1j for i in range(self.nmodes)] for j in range(self.nmodes)])
		for i in range(self.nmodes):
			for j in range(self.nmodes):
				rho[i][j] = common_factor * (0.5 * SIL[j][i] + FCE[i] * SIAP[j])
		return(rho, sqrt_updated)
	def integrand(self):
		rho_set = self.rho()
		rho = rho_set[0]
		sqrt_updated = rho_set[1]
		product = 0. + 0. * 1j
		for i in range(self.nmodes):
			for j in range(self.nmodes):
				product += self.NAC_matrix[i][j] * rho[i][j]
		integrand = np.exp(- 1j * self.adiabatic * self.time) * product
		return (integrand, sqrt_updated, rho)

correlation = []
trace_vector = []
updated_sqrt = - 1.0 + 0. * 1j
timestep = 0.01
for i in range(200000):
	time = timestep * i
	#if i // 1000 == 0:
	#	print(i)
	test = IC(nmodes, duschinsky, displacement, freq_initial, freq_final, time, beta, updated_sqrt, partition, adiabatic, R_matrix).integrand()
	updated_sqrt = test[1]
	if i == 0:
		correlation.append(np.real(test[0]) * 0.5)
	else:
		correlation.append(np.real(test[0]))
	trace_vector.append(np.trace(test[2]))

correlation = np.array(correlation)
trace_vector = np.array(trace_vector)
np.savetxt('integrand.txt', correlation)
np.savetxt('trace.txt', trace_vector)
print(sum(correlation) * timestep / atomic_time * 2)
