#################################################
# CCD (Couple-Cluster Doubles)
# Julie Butler Hartley
# Version 0.0.5
# Date Created: November 30, 2020
# Last Modified: November 30, 2020
#
# Python version of the couple-cluster doubles code algorithm presented
# in Chapter 8 of LNP 936
# Code assumes two-body matrix elements and single particle energies are
# given in lists
#
# TO-DOs: 
# Review C++ and make a C++ version of this file
# Find examples of TBME lists to test the code with
#################################################

#################################################
# OUTLINE
#################################################
# index(p, q, r, s, N_sp): Calculates and returns an index for t_amplitude
# and TBME lists that corresponds to the four single particle states passed
# in the arguments.
#
# CCD_Corr_Energy(t_amplitudes, N_Fermi, N_sp): Calculates the couple-cluster
# doubles correlation energy.
#
# t_amplitudes_zero(N_sp): 		Calculates the initial t amplitudes for a
# couple-cluster doubles iterative calculation.
#
# CCD_Update(t_amplitides, N_Fermi, N_sp)
#
# CCD_Iterate (TBME, SP_energy, N_Fermi, N_sp)

#############################
# IMPORTS
#############################
# NONE.

# INDEX
def index(p, q, r, s, N_sp):
	"""
		Inputs:
			p, q, r, s (ints): labels for single particle states
			N_sp (an int): the number of single particle states
		Returns:
			Unnamed (an int): The calculated index to extract information
				from lists for the specific four single particle states
		Calculates and returns an index for t_amplitude and TBME lists that
		corresponds to the four single particle states passed in the 
		arguments.		
	"""
	return p*N_sp**3 + q*N_sp**2 + r*N_sp + s

# CCD_CORR_ENERGY
def CCD_Corr_Energy(t_amplitudes, N_Fermi, N_sp):
	"""
		Inputs:
			t_amplitudes (a list): the t amplitudes from the couple cluster
				doubles calculation
			N_Fermi (an int): the number of single particle states below
				the Fermi level
			N_sp (an int): the number of single particle states
		Returns:
			energy (a double): the couple-cluster doubles correlation energy
		Calculates the couple-cluster doubles correlation energy.
	"""
	energy = 0
	# i and j are below the Fermi level
	for i in range (0, N_Fermi):
		for j in range (0, N_Fermi):
			# a and b are above the Fermi level
			for a in range (N_Fermi, N_sp):
				for b in range (N_Fermi, N_sp):
					# Perform the sum that calcualtes the correlation energy
					energy = energy + 0.25*TBME[index(i,j,a,b, N_sp)]*\
						t_amplitudes[index(a,b,i,j,N_sp)]
	# Return the correlation energy
	return energy

# T_AMPLITUDES_ZERO
def t_amplitudes_zero(N_sp):
	"""
		Inputs:
		Returns
		Calculates the initial t amplitudes for a couple-cluster doubles
		iterative calculation.
	"""
	# WHAT WOULD THIS VALUE BE????????????
	t_amplitudes = np.zeros(len(....))
	for i in range (0, N_Fermi):
		for j in range (0, N_Fermi):
			for a in range (N_Fermi, N_sp):
				for b in range (N_Fermi, N_sp):					
					t_amplitudes[index(a,b,i,j,N_sp)] = TBME[index(a,b,i,j,N_sp)]/(SP_energy[i] + SP_energy[j] - SP_energy[a] - SP_energy[b])
	return t_amplitudes	

# CCD_UPDATE
def CCD_Update(t_amplitides, N_Fermi, N_sp):
	t_amplitudes = np.zeros(len(t_amplitudes_old))
	# i and j are below the Fermi level
	for i in range (0, N_Fermi):
		for j in range (0, N_Fermi):
			# a, b, c, and d are above the Fermi level
			for a in range (N_Fermi, N_sp):
				for b in range (N_Fermi, N_sp):
					P_ijab = P(i,j,a,b)
					P_ij = P(i,j)
					P_ab = P(a,b)
					# Term 1
					# Calculate the two body matrix element
					sum = TBME[index(a, b, i, j, N_sp)]
					# Term 2
					for c in range (N_Fermi, N_sp):
						for d in range (N_Fermi, N_sp):
							# Add in contributions from particles c and d
							sum = sum + 0.5*TBME[index(a,b,c,d, N_sp)]*t_amplitudes_old[index(c,d,i,j, N_sp)]
					
					# Term 3
					for k in range (0, N_Fermi):
						for l in range (0, N_Fermi):
							sum = sum + 0.5*TBME[index(k,l,i,j,N_sp)]*t_amplitudes_old[index(a,b,k,l,N_sp)]

					# Term 4
					for k in range (0, N_Fermi):
						for c in range (N_Fermi, N_sp):
							sum = P_ijab*TBME[index(k,b,c,j, N_sp)]*t_amplitudes_old[index(a,c,i,k,N_sp)]

					# Terms 5-8
					for k in range (0, N_Fermi):
						for l in range (0, N_Fermi):
							for c in range (N_Fermi, N_sp):
								for d in range (N_Fermi, N_sp):
									# Term 5
									sum = sum + 0.25*TBME[index(k,l,c,d,N_sp)]*t_amplitudes_old[index(c,d,i,j,N_sp)]*t_amplitudes_old[index(a,b,k,l,N_sp)]
									# Term 6
									sum = sum + P_ij*TBME[index(k,l,c,d,N_sp)]*t_amplitudes_old[index(a,c,i,k,N_sp)]*t_amplitudes_old[index(b,d,j,l,N_sp)]
									# Term 7
									sum = sum - 0.5*P_ij*TBME[index(k,l,c,d,N_sp)]*t_amplitudes_old[index(d,c,i,k,N_sp)]*t_amplitudes_old[index(a,b,l,j,N_sp)]
									# Term 8
									sum = sum - 0.5*P_ab*TBME[index(k,l,c,d,N_sp)]*t_amplitudes_old[index(a,c,l,k,N_sp)]*t_amplitudes_old[index(d,b,i,j,N_sp)]
					
					# Calculate the energy denominator
					energy_demon = SP_energy[i] + SP_energy[j] - SP_energy[a] - SP_energy[b]
					
					# Find the new t amplitudes
					t_amplitudes[index(a,b,i,j)] = sum/energy_demon
	return t_amplitudes

# CCD_ITERATE
def CCD_Iterate (TBME, SP_energy, N_Fermi, N_sp):
	# Calculate the initial t amplitudes and correlation energy
	t_amplitudes_old = t_amplitudes_zero(N_sp)
	correlation_Energy_old = CCD_Corr_Energy(t_amplitudes_old, N_Fermi, N_sp)

	# Set the tolerance
	tolerance = 1e-5

	# Iterate until the correlation energy converges
	while (abs(energy_Diff) > tolerance):
		# Calculate the new t amplitudes
		t_amplitudes = CCD_Update(t_amplitides_old, N_Fermi, N_sp)
		# Calculate the new correlation energy
		correlation_Energy = CCD_Corr_Energy(t_amplitudes, N_Fermi, N_sp)
		# Calculate the energy difference and set up for next loop
		energy_Diff = correlation_Energy - correlation_Energy_old
		correlation_Energy_old = correlation_Energy
		t_amplitudes_old = t_amplitudes

	# Print converged correlation energy	
	print('Final Correlation Energy: ', correlation_Energy)				