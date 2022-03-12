import numpy as np
import datetime
import time
import pickle
import pdb
import gzip
import os
import sys
import sys
import re

#Declare the author that corresponds to the phenotype you want
#Also the indices of the group we want phenotypes for
author = sys.argv[2]
start_ind = int(sys.argv[1])
print(author)
print(start_ind)
#sys.exit()
#end_ind = int(sys.argv[2])
end_ind = start_ind + 50000

#Function to read in the data
def normRead(fileName, withHeader = True, delim = '\t', removeQuote = False):
	with open(fileName,"r",encoding = "latin-1") as f:
		totalData=[]
		for line in f.read().splitlines():
			if removeQuote:
				totalData.append(line.replace('"', '').strip().split(delim))
			else:
				totalData.append(line.split(delim))
	if withHeader:
		header=totalData[0]
		del totalData[0]
	else:
		header = None
	totalData=np.array(totalData)
	return(totalData,header)

#Read in the files that were split from the main UKBB phenotype file
date_ass,ass_header = normRead("date_assessed.csv", True, ",", True)
self_rep, self_rep_header = normRead("self_report_diag.csv", True, ",", True)
big_eid, eid_head = normRead("eid.csv", True, ",", True)

#sort the files just read in so they are in EID order
date_ass = date_ass[big_eid[:,0].argsort(),:]
self_rep = self_rep[big_eid[:,0].argsort(),:]
big_eid = big_eid[big_eid[:,0].argsort(),:]

#Get the phase, or the assessment type for cancer, self-report and medication report
#Phase meaning what index, or what time the person came in to be assessed (first, second or third assessment)
self_phase = [x.split("-")[1].split(".")[0] for x in self_rep_header]

#Read in the hesin type of data and again sort by EID
diag, head_diag = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt")
oper, head_oper = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_oper.txt")
hesin, head_hesin = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin.txt")
death, head_death = normRead("/home/kulmsc/athena/ukbiobank/hesin/death.txt")
cod, head_cod = normRead("/home/kulmsc/athena/ukbiobank/hesin/death_cause.txt")
diag = diag[diag[:,0].argsort(),:]
oper = oper[oper[:,0].argsort(),:]
hesin = hesin[hesin[:,0].argsort(),:]
death = death[death[:,0].argsort(),:]
cod = cod[cod[:,0].argsort(),:]
diag_eid = np.unique(diag[:,0])
oper_eid = np.unique(diag[:,0])
death_eid = np.unique(death[:,0])
cod_eid = np.unique(cod[:,0])
diag_max = diag.shape[0] - 1
oper_max = oper.shape[0] - 1
death_max = death.shape[0] - 1
cod_max = cod.shape[0] - 1
diag_eid_list = diag[:,0].tolist()
oper_eid_list = oper[:,0].tolist()
death_eid_list = death[:,0].tolist()
cod_eid_list = cod[:,0].tolist()

#Read in the disease defintion and split the types with multiple condictions, then make into a dict
defs, head_defs = normRead("fresh_jc_codes.txt")
def_ind = defs[:,0].tolist().index(author)
split_defs = list()
for type_def in defs[def_ind, 3:9]:
	split_defs.append(type_def.split("|"))

#Produce a dictionary so we can easily pull out the codes for each data type
use_defs = dict(zip(head_defs[3:9], split_defs))

#This if for timing, which I do not do anymore
#cancer_sr_time = 0
#noncancer_sr_time = 0
#hesin_setup_time = 0
#icd9_time = 0
#icd10_time = 0
#hesin_oper_time = 0
#med_time = 0


#Find the splits of the input data, so RAM will be less
#Simply find the eid corresponding to the start and stop index
#Continuing to skip them if they do not appear in the given data type
def subset_big_data(eid_list, big_data):

	for ind in range(start_ind, end_ind):
		if big_eid[ind][0] in eid_list:
			start_eid = eid_list.index(big_eid[ind][0])
			break

	for ind in range(end_ind, start_ind,-1):
		if big_eid[ind][0] in eid_list:
			end_eid = np.where(big_data[:,0] == big_eid[ind][0])[0][-1]
			break

	return(big_data[start_eid:(end_eid+1),:])
	

#I'm still bad at python indexing so I have to include this line to fix a single case example
if end_ind > big_eid.shape[0]:
	end_ind = big_eid.shape[0] - 1


#Apply the subset_big_data, and convert some things to lists for easier indexing
diag = subset_big_data(diag_eid_list, diag)
oper = subset_big_data(oper_eid_list, oper)
hesin = subset_big_data(hesin[:,0].tolist(), hesin)
diag_eid_list = diag[:,0].tolist()
oper_eid_list = oper[:,0].tolist()

#Establish indicdes of diag and oper for easy indexing (do not need to look and make index every time
#only need to add the extent of the current eid to this index)
start_eid_diag = 0
start_eid_oper = 0

#Subset the things that do not require the subset_big_data function
self_rep = self_rep[start_ind:end_ind,:]
date_ass = date_ass[start_ind:end_ind,:]
big_eid = big_eid[start_ind:end_ind,:]


#Prepare the data objects that will hold diagnosis results
df_occur = np.zeros((big_eid.shape[0], 5))
df_date = np.tile("__________", (big_eid.shape[0], 5))

#Iterate through the indices of each unique EID
#Note for lessons learned in speed-up:
# - do not subset objects (I think copies are made in RAM which really slow things down)
# - always try to just iterate to the next index rather than searching for some keyword
for ind in range(0, (end_ind - start_ind)):
	if ind % 10000 == 0:
		print(ind) #So I know how fast things are running

	#if ind > 29530:
	#	pdb.set_trace()

	#Set the value of the actual eid
	curr_eid = big_eid[ind][0]


	#NON-CANCER SELF REPORT #this is the same process as with the cancer data, but now with non-cancer data
	if use_defs["noncancer"][0] != "NA":
		#t22 = time.time()
		if any(np.isin(use_defs["noncancer"], self_rep[ind,:])):
			locs = np.where(np.isin(self_rep[ind,:], use_defs["noncancer"]))[0]
			try_date = date_ass[ind, int(self_phase[locs[0]])]
			if try_date == "":
				try_date = date_ass[ind, 0]
			df_date[ind,0] = try_date
			df_occur[ind,0] = 1

		#t3 = time.time()
		#noncancer_sr_time += t3-t22


	#HESIN DIAG
	if curr_eid in diag_eid_list: #check if the eid has ever had any ICD codes
		#t4 = time.time()
		#set up the hesin indexing, since multiple rows can have the same eid we use this
		#use this technique to get the index of the next unique EID
		if curr_eid != diag_eid_list[-1]:
			next_eid = np.unique(diag[start_eid_diag:(start_eid_diag+7000),0])[1]
			ext_eid = diag_eid_list.index(next_eid)
		else:
			ext_eid = diag_max

		#t5 = time.time()
		#hesin_setup_time += t5-t4

		#ICD - 9
		#as was explained either the full or slightly shorter ICD code may match what we want, so we create both
		#use_short_icd9_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,4]]
		#Then we check to see if either the short short or long ICD codes (from the hesin) match any of the definitions
		#if any(np.isin(use_defs["icd9"], use_short_icd9_diag)):
		#if any(np.isin(use_defs["icd9"], diag[start_eid_diag:ext_eid,4])):
		#	#If there is a match we get all of the locations in the hesin of the match
		#	short_icd9_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,4], use_defs["icd9"]))[0][0]
		#else:
		short_icd9_locs = 10000
		if any(np.isin(use_defs["icd9"], diag[start_eid_diag:ext_eid,4])):
                        long_icd9_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,4], use_defs["icd9"]))[0][0]
		else:
			long_icd9_locs = 10000
		#There cannot possibly be a match greater than 10000, so I set it to the impossible value
		if long_icd9_locs != 10000 or short_icd9_locs != 10000:
			#Find the minimum matching location as it is the earliest time
			icd9_locs = min((long_icd9_locs, short_icd9_locs))
			#Then we get the instance or the numbered time the EID went to the hospital
			icd9_ins_index = diag[start_eid_diag+icd9_locs,1]

			#We carry the ins_index and EID to the larger hesin array to get the date
			#Somtimes the most accurate form of the date is empty so we go to the next best
			raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),5][0]
			if raw_date == "":
				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),21][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),4][0]

			df_occur[ind,1] = 1
			df_date[ind,1] = raw_date
		
		#t6 = time.time()
		#icd9_time += t6-t5

		#The process for ICD10 is nearly exactly the same as for ICD9
		#ICD - 10
		use_short_icd10_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,6]]
		if any(np.isin(use_defs["icd10"], use_short_icd10_diag)):
			short_icd10_locs = np.where(np.isin(use_short_icd10_diag, use_defs["icd10"]))[0][0]
		#if any(np.isin(use_defs["icd10"], diag[start_eid_diag:ext_eid,6])):
		#	short_icd10_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,6], use_defs["icd10"]))[0][0]
		else:
			short_icd10_locs = 10000
		if any(np.isin(use_defs["icd10"], diag[start_eid_diag:ext_eid,6])):
			long_icd10_locs = np.where(np.isin(diag[start_eid_diag:ext_eid,6], use_defs["icd10"]))[0][0]
		else:
			long_icd10_locs = 10000
		if long_icd10_locs != 10000 or short_icd10_locs != 10000:
			icd10_locs = min((long_icd10_locs, short_icd10_locs))
			icd10_ins_index = diag[start_eid_diag+icd10_locs,1]

			raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),5][0]
			if raw_date == "":
				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),21][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),4][0]

			df_occur[ind,2] = 1
			df_date[ind,2] = raw_date
		start_eid_diag = ext_eid

		#t7 = time.time()
		#icd10_time += t7-t6


	#The process for HESIN OPER is nearly exactly the same as for HESIN DIAG (which includes ICD10 and ICD))
	#HESIN OPER
	if curr_eid in oper_eid_list and use_defs["opcs"][0] != "NA":
		#t8 = time.time()
		if curr_eid != oper_eid_list[-1]: #if it is not the last in the list
			next_eid = np.unique(oper[start_eid_oper:(start_eid_oper+7000),0])[1]
			ext_eid = oper_eid_list.index(next_eid)
		else:
			ext_eid = oper_max

		#use_short_oper = [x[0:3] for x in oper[start_eid_oper:ext_eid,7]]
		#if any(np.isin(use_defs["opcs"], use_short_oper)):
		if any(np.isin(use_defs["opcs"], oper[start_eid_oper:ext_eid,7])):
			short_oper_locs = np.where(np.isin(oper[start_eid_oper:ext_eid,7], use_defs["opcs"]))[0][0]
		else:
			short_oper_locs = 10000
		if any(np.isin(use_defs["opcs"], oper[start_eid_oper:ext_eid,7])):
			long_oper_locs = np.where(np.isin(oper[start_eid_oper:ext_eid,7], use_defs["opcs"]))[0][0]
		else:
			long_oper_locs = 10000
		if long_oper_locs != 10000 or short_oper_locs != 10000:
			oper_locs = min((long_oper_locs, short_oper_locs))
			oper_ins_index = oper[start_eid_oper+oper_locs,1]

			raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),5][0]
			if raw_date == "":
				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),21][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),4][0]

			df_occur[ind,3] = 1
			df_date[ind,3] = raw_date
		start_eid_oper = ext_eid

		#t9 = time.time()
		#hesin_oper_time += t9-t8

	#The medication process is nearly the same as for cancer and non-cancer
	#MEDICATION
	if use_defs["icd10"][0] != "NA":
		if curr_eid in death_eid_list:
			if np.any(np.isin(cod[cod[:,0] == curr_eid,4], use_defs["icd10"])):
				df_occur[ind,4] = 1
				df_date[ind,4] = death[death_eid_list.index(curr_eid),4]
				



#save the output phenos
np.savetxt("raw_output/diag.coding." + author.lower() + "." + str(start_ind) + ".txt.gz", df_occur, fmt='%i')
np.savetxt("raw_output/diag.time." + author.lower() + "." + str(start_ind) + ".txt.gz", df_date, fmt='%s')



print("end")
