
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
from operator import itemgetter

def plot_results(target, prediction, filename, R2 = 666):
	plt.plot(target, prediction, 'ro')
	x=[max([min(target)-0.01, min(prediction)-0.01]), min([max(target)+0.01,max(prediction)+0.01])]
	y=x
	plt.plot(x, y, 'k:', lw=2)
	plt.xlabel("Target")
	plt.ylabel("Predictions")
	plt.title("Target vs Prediction: R2 = " + str(round(R2,4)))
	plt.savefig(filename)

def plot_controls(dataf, filename, n_columns=5):
	trans=dataf.T
	pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
	for c in trans.columns:
		control = dict(trans[c])
		control_biggest = dict(sorted(control.items(), key = itemgetter(1), reverse = True)[:n_columns])
		fig,ax = plt.subplots()  #create a new figure
		ax.bar(list(control_biggest.keys()), list(control_biggest.values()))
		ax.set_title(c)
		ax.tick_params(labelrotation=12, labelsize=10)
		fig.savefig(pdf, format='pdf')
	pdf.close()
	return(control_biggest)