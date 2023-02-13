
import matplotlib.pyplot as plt

def plot_results(target, prediction, filename):
	plt.plot(target, prediction, 'ro')
	x=[max([min(target)-0.01, min(prediction)-0.01]), min([max(target)+0.01,max(prediction)+0.01])]
	y=x
	plt.plot(x, y, 'k:', lw=2)
	plt.savefig(filename)
