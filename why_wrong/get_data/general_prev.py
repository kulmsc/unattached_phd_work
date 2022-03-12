import numpy as np
from datetime import datetime
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
author = "Bentham" 
start_ind = int(sys.argv[1])
end_ind = start_ind + 50000 #int(sys.argv[2])

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
date_ass,ass_header = normRead("/home/kulmsc/athena/doc_score/analyze_score/construct_defs/date_assessed.csv", True, ",", True)
meds, meds_header = normRead("/home/kulmsc/athena/doc_score/analyze_score/construct_defs/medications.csv", True, ",", True)
self_rep, self_rep_header = normRead("/home/kulmsc/athena/doc_score/analyze_score/construct_defs/self_report_diag.csv", True, ",", True)
cancer_rep, cancer_rep_header = normRead("/home/kulmsc/athena/doc_score/analyze_score/construct_defs/self_report_cancer.csv", True, ",", True)
big_eid, eid_head = normRead("/home/kulmsc/athena/doc_score/analyze_score/construct_defs/eid.csv", True, ",", True)

#sort the files just read in so they are in EID order
date_ass = date_ass[big_eid[:,0].argsort(),:]
meds = meds[big_eid[:,0].argsort(),:]
self_rep = self_rep[big_eid[:,0].argsort(),:]
cancer_rep = cancer_rep[big_eid[:,0].argsort(),:]
big_eid = big_eid[big_eid[:,0].argsort(),:]
np.savetxt("big_eid.txt", big_eid, fmt='%s')


#Get the phase, or the assessment type for cancer, self-report and medication report
#Phase meaning what index, or what time the person came in to be assessed (first, second or third assessment)
cancer_phase = [x.split("-")[1].split(".")[0] for x in cancer_rep_header]
self_phase = [x.split("-")[1].split(".")[0] for x in self_rep_header]
meds_phase = [x.split("-")[1].split(".")[0] for x in meds_header]

#Read in the hesin type of data and again sort by EID
diag, head_diag = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt")
oper, head_oper = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin_oper.txt")
hesin, head_hesin = normRead("/home/kulmsc/athena/ukbiobank/hesin/hesin.txt")
diag = diag[diag[:,0].argsort(),:]
oper = oper[oper[:,0].argsort(),:]
hesin = hesin[hesin[:,0].argsort(),:]
diag_eid = np.unique(diag[:,0])
oper_eid = np.unique(diag[:,0])
diag_max = diag.shape[0] - 1
oper_max = oper.shape[0] - 1
diag_eid_list = diag[:,0].tolist()
oper_eid_list = oper[:,0].tolist()

poss_diag_icd9 = np.unique(np.array([x[0:3] for x in diag[:,4]]))
poss_diag_icd10 = np.unique(np.array([x[0:3] for x in diag[:,6]]))
poss_oper = np.unique(np.array([x[0:3] for x in oper[:,7]]))
poss_meds = np.unique(meds)
poss_cancer = np.unique(cancer_rep)
poss_self = np.unique(self_rep)

poss_diag_icd9 = poss_diag_icd9[poss_diag_icd9 != ""]
poss_diag_icd10 = poss_diag_icd10[poss_diag_icd10 != ""]
poss_oper = poss_oper[poss_oper != ""]
poss_meds = poss_meds[poss_meds != ""]
poss_cancer = poss_cancer[poss_cancer != ""]
poss_self = poss_self[poss_self != ""]

#Read in the disease defintion and split the types with multiple condictions, then make into a dict
defs, head_defs = normRead("/home/kulmsc/athena/doc_score/analyze_score/descript_defs/author_defs")
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
meds = meds[start_ind:end_ind,:]
self_rep = self_rep[start_ind:end_ind,:]
date_ass = date_ass[start_ind:end_ind,:]
cancer_rep = cancer_rep[start_ind:end_ind,:]
big_eid = big_eid[start_ind:end_ind,:]

real_use_date = [datetime.strptime(x, "%Y-%m-%d") for x in date_ass[:,0]]

#Prepare the data objects that will hold diagnosis results
df_occur = np.zeros((big_eid.shape[0], 6))
df_date = np.tile("__________", (big_eid.shape[0], 6))
df_ass_status = np.zeros((big_eid.shape[0], 6))

df_prev_icd10 =  np.zeros((big_eid.shape[0], len(poss_diag_icd10)))
df_prev_icd9 = np.zeros((big_eid.shape[0], len(poss_diag_icd9)))
df_prev_cancer = np.zeros((big_eid.shape[0], len(poss_cancer)))
df_prev_noncancer = np.zeros((big_eid.shape[0], len(poss_self)))
df_prev_opcs = np.zeros((big_eid.shape[0], len(poss_oper)))
df_prev_meds = np.zeros((big_eid.shape[0], len(poss_meds)))


df_inci_icd10 =  np.zeros((big_eid.shape[0], len(poss_diag_icd10)))
df_inci_icd9 = np.zeros((big_eid.shape[0], len(poss_diag_icd9)))
df_inci_opcs = np.zeros((big_eid.shape[0], len(poss_oper)))

df_date_icd10 =  np.tile(' ', (big_eid.shape[0], len(poss_diag_icd10)))
df_date_icd9 = np.tile(' ', (big_eid.shape[0], len(poss_diag_icd9)))
df_date_opcs = np.tile(' ', (big_eid.shape[0], len(poss_oper)))

cancer_rep = cancer_rep[:,cancer_phase == '0']
self_rep = self_rep[:,self_phase == '0']

#Iterate through the indices of each unique EID
#Note for lessons learned in speed-up:
# - do not subset objects (I think copies are made in RAM which really slow things down)
# - always try to just iterate to the next index rather than searching for some keyword
for ind in range(0, (end_ind - start_ind)):
	#print(ind)
	if ind % 10000 == 0:
		print(ind) #So I know how fast things are running


	#Set the value of the actual eid
	curr_eid = big_eid[ind][0]

	#Now we go through and fill in each column of df_occur and df_date, one data type at a time

	#CANCER SELF REPORT (this is the name of the column in df_occur and df_date I am filling in)
	if any(cancer_rep[ind,:] != ""): #check to see if there is a noncancer definition
		#t1 = time.time()
		for cancer_code in cancer_rep[ind,cancer_rep[ind,:] != ""]:
			df_prev_cancer[ind, poss_cancer == cancer_code] = 1

		#t2 = time.time()
		#cancer_sr_time += t2-t1

	#NON-CANCER SELF REPORT #this is the same process as with the cancer data, but now with non-cancer data
	if any(self_rep[ind,:] != ""):
		#t22 = time.time()
		for self_code in self_rep[ind, self_rep[ind,:] != ""]:
			df_prev_noncancer[ind, poss_self == self_rep] = 1
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


		use_date = date_ass[ind,0]


		#ICD - 9
		#as was explained either the full or slightly shorter ICD code may match what we want, so we create both
		#use_short_icd9_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,4]]
		if np.any(diag[start_eid_diag:ext_eid,4] != ""):
			use_short_icd9_diag = np.array([x[0:3] for x in diag[start_eid_diag:ext_eid,4]])
			short_icd9_locs = np.where(np.invert(np.in1d(use_short_icd9_diag, "")))[0]

			for kk in range(len(short_icd9_locs)):

				#Then we get the instance or the numbered time the EID went to the hospital
				icd9_ins_index = diag[start_eid_diag+short_icd9_locs[kk],1]

				#We carry the ins_index and EID to the larger hesin array to get the date
				#Somtimes the most accurate form of the date is empty so we go to the next best
				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),5][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),21][0]
					if raw_date == "":
						raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd9_ins_index),4][0]
						if raw_date == "":
							raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")
				try:
					raw_date = datetime.strptime(raw_date, "%d/%m/%Y")
				except:
					raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")
				if raw_date < real_use_date[ind]:
					df_prev_icd9[ind, poss_diag_icd9 == use_short_icd9_diag[short_icd9_locs[kk]]] = 1
				df_inci_icd9[ind, poss_diag_icd9 == use_short_icd9_diag[short_icd9_locs[kk]]] = 1
				df_date_icd9[ind, poss_diag_icd9 == use_short_icd9_diag[short_icd9_locs[kk]]] = str(raw_date).split(' ')[0]


		
		#t6 = time.time()
		#icd9_time += t6-t5

		#The process for ICD10 is nearly exactly the same as for ICD9
		#ICD - 10
		#use_short_icd10_diag = [x[0:3] for x in diag[start_eid_diag:ext_eid,6]]
		if np.any(diag[start_eid_diag:ext_eid,6] != ""):
			use_short_icd10_diag = np.array([x[0:3] for x in diag[start_eid_diag:ext_eid,6]])
			short_icd10_locs = np.where(np.invert(np.in1d(use_short_icd10_diag, "")))[0]

			for kk in range(len(short_icd10_locs)):

                                #Then we get the instance or the numbered time the EID went to the hospital
				icd10_ins_index = diag[start_eid_diag+short_icd10_locs[kk],1]

				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),5][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),21][0]
					if raw_date == "":
						raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == icd10_ins_index),4][0]
						if raw_date == "":
							raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")


				try:
					raw_date = datetime.strptime(raw_date, "%d/%m/%Y")
				except:
					raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")
				if raw_date < real_use_date[ind]:
					df_prev_icd10[ind, poss_diag_icd10 == use_short_icd10_diag[short_icd10_locs[kk]]] = 1
				df_date_icd10[ind, poss_diag_icd10 == use_short_icd10_diag[short_icd10_locs[kk]]] = str(raw_date).split(' ')[0]
				df_inci_icd10[ind, poss_diag_icd10 == use_short_icd10_diag[short_icd10_locs[kk]]] = 1


		start_eid_diag = ext_eid

		#t7 = time.time()
		#icd10_time += t7-t6


	#The process for HESIN OPER is nearly exactly the same as for HESIN DIAG (which includes ICD10 and ICD))
	#HESIN OPER

	if curr_eid in oper_eid_list:
		#t8 = time.time()
		if curr_eid != oper_eid_list[-1]: #if it is not the last in the list
			next_eid = np.unique(oper[start_eid_oper:(start_eid_oper+7000),0])[1]
			ext_eid = oper_eid_list.index(next_eid)
		else:
			ext_eid = oper_max

		#use_short_oper = [x[0:3] for x in oper[start_eid_oper:ext_eid,7]]
		if np.any(oper[start_eid_oper:ext_eid,7] != ""):
			use_short_oper = np.array([x[0:3] for x in oper[start_eid_oper:ext_eid,7]])
			short_oper_locs = np.where(use_short_oper != "")[0]

			for kk in range(len(short_oper_locs)):

				oper_ins_index = oper[start_eid_oper+short_oper_locs[kk],1]

				raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),5][0]
				if raw_date == "":
					raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),21][0]
					if raw_date == "":
						raw_date = hesin[np.logical_and(hesin[:,0] == curr_eid, hesin[:,1] == oper_ins_index),4][0]
						if raw_date == "":
							raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")
				
				try:
					raw_date = datetime.strptime(raw_date, "%d/%m/%Y")
				except:
					raw_date = datetime.strptime("01/01/2021", "%d/%m/%Y")

				#if curr_eid == "2505970":
				#	pdb.set_trace()

				if raw_date < real_use_date[ind]:
					#if poss_oper[poss_oper == use_short_oper[kk]] == "B27":
					#	pdb.set_trace()
					#if poss_oper[poss_oper == use_short_oper[kk]] == "B28":
					#	pdb.set_trace()

					df_prev_opcs[ind, poss_oper == use_short_oper[short_oper_locs[kk]]] = 1
				df_inci_opcs[ind, poss_oper == use_short_oper[short_oper_locs[kk]]] = 1
				df_date_opcs[ind, poss_oper == use_short_oper[short_oper_locs[kk]]] = str(raw_date).split(' ')[0]


		start_eid_oper = ext_eid

		#t9 = time.time()
		#hesin_oper_time += t9-t8

	#The medication process is nearly the same as for cancer and non-cancer
	#MEDICATION
	#t10 = time.time()
	if any(meds[ind,:] != ""):

		for med_code in meds[ind ,meds[ind,:] != ""]:
                        df_prev_meds[ind, poss_meds == med_code] = 1

		#t11 = time.time()
		#med_time += t11 - t10



#save the output phenos
np.savetxt("raw_output/diag.icd10." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_icd10, fmt='%i')
np.savetxt("raw_output/diag.icd9." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_icd9, fmt='%i')
np.savetxt("raw_output/diag.cancer." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_cancer, fmt='%i')
np.savetxt("raw_output/diag.noncancer." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_noncancer, fmt='%i')
np.savetxt("raw_output/diag.opcs." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_opcs, fmt='%i')
np.savetxt("raw_output/diag.meds." + author.lower() + "." + str(start_ind) + ".txt.gz", df_prev_meds, fmt='%i')


np.savetxt("raw_output/diag.icd10." + author.lower() + ".txt.gz", poss_diag_icd10, fmt='%s')
np.savetxt("raw_output/diag.icd9." + author.lower() + ".txt.gz", poss_diag_icd9, fmt='%s')
np.savetxt("raw_output/diag.cancer." + author.lower() + ".txt.gz", poss_cancer, fmt='%s')
np.savetxt("raw_output/diag.noncancer." + author.lower() + ".txt.gz", poss_self, fmt='%s')
np.savetxt("raw_output/diag.opcs." + author.lower() + ".txt.gz", poss_oper, fmt='%s')
np.savetxt("raw_output/diag.meds." + author.lower() + ".txt.gz", poss_meds, fmt='%s')


np.savetxt("raw_output/inci.icd10." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_icd10, fmt='%i')
np.savetxt("raw_output/inci.icd9." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_icd9, fmt='%i')
np.savetxt("raw_output/inci.opcs." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_opcs, fmt='%i')

np.savetxt("raw_output/date.icd10." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_icd10, fmt='%s')
np.savetxt("raw_output/date.icd9." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_icd9, fmt='%s')
np.savetxt("raw_output/date.opcs." + author.lower() + "." + str(start_ind) + ".txt.gz", df_inci_opcs, fmt='%s')


print("end")
