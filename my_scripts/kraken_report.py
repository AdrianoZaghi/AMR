import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="Process some integers.")
parser.add_argument('filename', type=str, help='The filename to process')
parser.add_argument('-o', '--output', type=str, help='The output filename')
args = parser.parse_args()


binning = {}
with open(args.filename, 'r') as f:
	lines = f.readlines()
	for line in lines:
		bin = line.split('\t')[2]
		if bin in binning.keys():
			binning[bin] += 1
		else:
			binning[bin] = 1

print("Report completato")

binning = pd.DataFrame.from_dict(list(binning.items()))
binning.columns = ["id", "count"]
with open(args.output, 'w') as f:
	binning.to_csv(f, sep='\t', header=False)
	print("\n\n", file=f)
	for i in range (0, 5):
		slice = binning[binning['id'].astype(str).str.startswith(str(i))]
		if slice.empty:
			print(str(i) + "   0", file=f)
		else:
			print(str(i) + "   " + str(slice["count"].sum()), file=f)

	