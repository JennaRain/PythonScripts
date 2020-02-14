# TE Landscapes
# Jenna R. Grimshaw
# August 11, 2019
# jenna.grimshaw81@gmail.com


# Arguments: input rm.bed file
# Arguments: How we determine what we are interested input
# Arguments: Output folder - Do we want to create a new directory?

#############################################
############## Import modules ###############
#############################################
import logging
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

##############################################
################## Arguments #################
##############################################

ap = argparse.ArgumentParser()
ap.add_argument("-b", "--BEDFILE", required=True, help="Please enter your bedfile")
ap.add_argument("-g", "--GENOMESIZE", required=True, type = int, help="Enter size of genome (bp)")
ap.add_argument("-t", "--TAXON", required=True, help="Enter the taxon abbreviation")
args=vars(ap.parse_args())

##############################################
################## Functions #################
##############################################

def SUMofSIZE(data, element, min, max):
	size = []
	for DIV, SIZE, TE in zip(data.Divergence, data.Size, data.TE):
		if element == TE:
			if DIV >= min and DIV < max:
				size.append(SIZE)
	return sum(size)	
# This function will add up the lengths/sizes 
# of TEs within a min and max Divergence
# and then will return the sum of those lengths



################################################
################################################
################################################

# Enter genome size
genomesize = args["GENOMESIZE"]
print("Genome size is: ", genomesize)

# Import Rm_bed file
#fields = ["Scaffold", "Start", "Stop", "TE", "Score", "Orientation", "Class", "Family", "Divergence"]
fields = ["Scaffold", "Start", "Stop", "TE", "Score", "Orientation", "Class", "Family", "Divergence", "XX"]
file = pd.read_csv(args["BEDFILE"], sep="\t", names=fields)
print(file.head())
print("##################")

# We are only interested in SINEs
# So create a new dataframe only using SINEs
# Calculate length, drop unneeded columns
# And reorder the columns
SINE = file[file.Class =="SINE"]
SINE["Size"] = SINE["Stop"] - SINE["Start"]
SINE=SINE.drop(["Start", "Stop", "Score", "Class", "Orientation", "Family", "XX"], axis=1)
SINE = SINE[["Scaffold", "Size", "TE","Divergence"]]
SINE["Divergence"] = SINE["Divergence"]/100	# If bedfile divergences are large numbers and need to be changed to decimals
print(SINE.head())
print("######################")

# Make sure everything is at least 100 bp long
print("minimum TE length", SINE["Size"].min()) 
print("######################")

# If you want to know which TE subfamilies have how many occurrences
#count=SINE.groupby("TE").count()
#print(count)
#count.to_csv(args["TAXON"]+"counts.csv")

# Keep only elements that occur at least 100 times
SINE = (SINE.groupby("TE").filter(lambda x : len(x)>100))
TELIST = (SINE["TE"].unique())
print("Filtered TE list: ", len(TELIST), TELIST)
print("#######################")

# Create dataframe to be used in loopfile
# I want one column with ranges of divergences
# Starting at 0 to 0.50 by increments of 0.01
loopfile = pd.DataFrame(columns = TELIST)
loopfile["series"] = range(51)
loopfile["Range"] = loopfile["series"]/100
loopfile = loopfile.drop("series", 1)
print("creating loopfile")
#print(loopfile) #empty loopfile to be filled in at next step
print("#######################")

# Run the program
for TE in TELIST:	# for every family
	min = 0					# Start at 0 min divergence
	data = []				# Empty list to be used in loop
	for element in loopfile.Range:	# For each value in Range
		if element == min/100:		# if value matches min/100
			#print("The min is: ",min)
			start = min/100			# Make min/100 the start divergence
			sizesum = SUMofSIZE(SINE, TE, start, start+0.01) # run function
			data.append(sizesum/genomesize) # Append the sum/genome size to get proportion
			min += 1		# Increase minimum by 1
	#print(data)
	loopfile[TE] = data # For that family, append the data to loopfile
print(loopfile)
loopfile.to_csv(args["TAXON"]+"-divergences.csv")
