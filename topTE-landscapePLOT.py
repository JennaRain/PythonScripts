# Check WIKI for example of appropriate input file

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-d", "--DIVFILE", required=True, help="a csv file that has sums for divergences 0.00-0.50")
ap.add_argument("-c", "--CUTOFF", default=15, type=int)
ap.add_argument("-t", "--TAXON", required=True, help="Taxon abbreviation to name figure")
args=vars(ap.parse_args())	

file = pd.read_csv(args["DIVFILE"], index_col=0)
print(file.head())

def gettop10(data, cutoff=args["CUTOFF"]):
	sums = data.sum(axis=0)
	sums = pd.DataFrame(sums)
	sums = sums.transpose().drop(["Range"], axis=1).sort_values(axis=1, by=0, ascending=False)
	print("**********************")
	print("sums", sums)
	tops = sums.iloc[:, 0:cutoff]
	top10 = list(tops)
	print(len(top10))
	print("**********************")
	return(top10)
#print(gettop10(file))


def keeptop10data(data, list):
	top10list = gettop10(file)
	top10data = pd.DataFrame()
	bottom = pd.DataFrame() 
	for element in file.columns:
		#print(element)
		if element in top10list:
			top10data[[element]] = file[[element]]
		else:
			bottom[[element]] = file[[element]]
	bottom = bottom.drop("Range", axis=1)
	print("Number of subfamilies", len(bottom.columns))
	bottom["Other"] = bottom.sum(axis=1)	
	top10data["Other (x)"] = bottom["Other"] #Substitute "x" for the number of subfamilies in "Other"
	top10data["Range"] = file["Range"]
	return(top10data)


top10elements = keeptop10data(file, gettop10(file,args["CUTOFF"]+1))
print(top10elements)
TELIST = list(top10elements)
TELIST = TELIST[0:args["CUTOFF"]+1]
print("Number of TEs", len(TELIST))
for TE in TELIST:
	print(TE)
	plt.plot("Range", TE, data = top10elements )
axes = plt.gca()
#axes.set_ylim([0,0.0023])	# Use if you want to set your y-axes to all match
plt.gca().get_lines()[7].set_color("darkblue")
plt.gca().get_lines()[10].set_color("sienna")
plt.gca().get_lines()[11].set_color("fuchsia")
plt.gca().get_lines()[12].set_color("crimson")
plt.gca().get_lines()[13].set_color("darkgreen")
plt.gca().get_lines()[14].set_color("plum")
plt.gca().get_lines()[-1].set_color("grey")
plt.ylabel("Proportion of Genome")
plt.xlabel("Divergence")
plt.legend()
plt.savefig("top"+args["TAXON"]+".png")
plt.show()
