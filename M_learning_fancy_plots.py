import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

cosa = "depth" # "design" or "depth"
if cosa == "design":
	R2s_decesare = pd.read_csv("M_learning_images/R2s_design_effects_decesare.csv")
	R2s_patric = pd.read_csv("M_learning_images/R2s_design_effects_patric.csv")

	fig, ax = plt.subplots()
	ax.errorbar(R2s_decesare['Number of species'], R2s_decesare['Average R2'], yerr=R2s_decesare['R2 errors'], fmt='o', color = 'r', label = "Assembled genomes")
	ax.errorbar(R2s_patric['Number of species'], R2s_patric['Average R2'], yerr=R2s_patric['R2 errors'], fmt='o', color = 'b', label = "PATRIC genomes")
	ax.set_xlabel('Number of species in the sample')
	ax.set_xticks([2,3,4,5])
	ax.axhline(y=1, color='r', linestyle='--')
	ax.set_ylabel('R2')
	ax.set_title("Effect of sample complexity")
	ax.set_ylim(0.5,1.1)
	ax.legend()

elif cosa == "depth":
	R2s_decesare = pd.read_csv("M_learning_images/R2s_depth_effects_decesare.csv")
	R2s_patric = pd.read_csv("M_learning_images/R2s_depth_effects_patric.csv")
	
	fig, ax = plt.subplots()
	ax.errorbar(R2s_decesare['Subsample'], R2s_decesare['Average R2'], yerr=R2s_decesare['R2 errors'], fmt='o', color = 'r', label = "Assembled genomes")
	ax.errorbar(R2s_patric['Subsample'], R2s_patric['Average R2'], yerr=R2s_patric['R2 errors'], fmt='o', color = 'b', label = "PATRIC genomes")
	plt.xscale('log')
	ax = plt.gca()
	ax.axhline(y=1, color='r', linestyle='--')
	ax.set_ylim(0,1.2)
	ax.set_title("Effect of sample depth")
	ax.set_xlabel('Sample size (Gbp)')
	ax.set_ylabel('R2')
	ax.legend()
