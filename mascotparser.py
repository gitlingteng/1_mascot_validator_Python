#mascotparser is a small program which parse the proteiomic mascot file line by line, save the peptide info
#dictionary object for return.
#


#!/usr/bin/env python
# coding: UTF-8
# enable debugging
import cgi, cgitb
import sys
import cPickle
import time
#import psycopg2
import bisect
import os
import itertools # built in to 2.6

import copy
import datetime
import time
import operator
from operator import itemgetter, attrgetter
import math
from random import shuffle, choice
#saved in local directory
import corestats # http://www.goldb.org/corestats.html
#from padnums import pprint_table
#sudo apt-get
from numpy import average
from numpy import std as stdev
from numpy import var as variance
from xlwt import Workbook,Style,easyxf # get these from http://www.python-excel.org/
import matplotlib.pyplot as pyplot
import matplotlib

global screen_print
global screen_print2
global screen_print3
global pickle_data
global dat_file_version
global proton_mass

proton_mass = 1.00727646688
dat_file_version ="Old" 

screen_print = True
screen_print3 =True
screen_print2 =True
single_run =True
pickled_data = False
#testre =test()


#main function, return a dictionary where the info of mascot file was saved
def test():    
   # filename = "/home/lteng/Desktop/Projects/6MassSpec/ling_code/databases/Mascot_file_Mann_HelaEpo_05_F2893.dat"
    filename ="Mascot_sim.dat"   #simple mascot dat file for test running
    pri_dic = getmassdict(filename)
    return pri_dic

#main sub function#
def getmassdict(filename):
        # print "Starting",filename
                current_dir = '/home/lteng/Desktop/Projects/6MassSpec/ling_code/'
                os.chdir(current_dir)
            #    filename="Mascot_sim.dat"
                data = readfile(filename)
                f = filename.split('/')[-1]
        #	configure output directory
                # os.mkdir(current_dir + '/'  + "single_output")
                # os.chdir(current_dir + '/'  + "single_output")
                # if not single_run and not os.path.exists(current_dir + '/' + 'output_files/'+data_folder + '/' + folder + '/' + f.split('.')[0]):
                # 	os.mkdir(current_dir + '/'  + 'output_files/'+data_folder + '/' + folder + '/' + f.split('.')[0])
                if not single_run:
                        os.chdir(current_dir + '/'  + 'output_files/'+data_folder + '/' + folder + '/' + f.split('.')[0])
                else:
                        os.chdir(current_dir + "/" + "single_output")
                new_data = []
                # load in the whole .DAT file
                for line in data:
                        if '\r' in line:
                                line = line[:-1]
                        new_data.append(line)
                data = new_data
                if screen_print:
                        print "Data file has",len(data),"lines."
                content = make_data_dictionary(data)
                sections = content.keys()	
                mass_data = content['masses']
                #print mass_data
                modifications = {0:['N',0.0]} # instantiate the dictionary of modifications
                # find the modification information in the DAT file
                for line in mass_data:
                        if 'delta' in line:
                                n = int(line.split('=')[0].split('delta')[1])
                                m = float(line.split('=')[1].split(',')[0])
                                aa = line.split(' ')[1].split('(')[1].split(')')[0]
                                if aa == "N-term":
                                        aa = '1'
                                elif aa == "K":
                                        lys_mod = n
                                        delta_lys6 = m
                                elif aa == "R":
                                        arg_mod = n
                                        delta_arg10 = m
                                modifications[n] = [aa,m]
                        '''attempt at adding fixed mods - this adds another dictionary item with the fixed mods under the key of "fixed" '''
                        # elif 'FixedMod' in line and 'Residues' not in line:
                        # 	aa = line.split(' ')[1].split('(')[1].split(')')[0]
                        # 	m = float(line.split('=')[1].split(',')[0])
                        # 	n = 'F' + aa
                        # 	modifications.setdefault('fixed',[]).append([aa,m])

                print modifications
                # will need to generalize this part for other isotope labeling
                # this script is currently only written for LYS/ARG SILAC
                # if not('lys_mod' in dir() or 'arg_mod' in dir()):
                # 	print filename, "not SILAC"
                # 	continue # with next file in list
                mods_to_check = [lys_mod, arg_mod]
                poss_mods = [delta_lys6, delta_arg10]
                MW_test_group = []
                mod_names = [["LYS",'K'],["ARG","R"]]
                mod_dict = {} # will print corresponding number for amino acid
                for x,m in enumerate([a[1] for a in mod_names]):
                        mod_dict[m] = mods_to_check[x]
                for m in "ABCDEFGHIJLMNOPQSTUVWXYZ":
                        mod_dict[m] = 0
                # possible combinations of mods allowed (e.g. ARG-LYS-LYS)
                combos = ['0','1','00','11','01']#,'011','110','000','111', '0111', '1000', '1100', '1111', '0000']
                for c in combos:
                        MW_test_group.append(["_".join([mod_names[int(b)][0] for b in c]),sum([poss_mods[int(cc)] for cc in c])])
                isotope_masses = {}
                # pre-calculate all the possible isotope masses
                for line in MW_test_group:
                        key = line[0]
                        value = line[1]
                        isotope_masses[key] = value
                queries = [content[a] for a in content if 'query' in a]
                # for q in queries[0:100]:
                # 	print q
                # exit()
                if not pickled_data: # if we are not using a smaller, pickled testing data set
                        query_dictionary = make_query_dictionary(content["summary"],content["peptides"],queries)
                        # for q in query_dictionary:
                        # 	print q
                        # 	print query_dictionary[q]
                        # 	exit()
                        if screen_print:
                                print "...done making target dictionary"			
                        keys = query_dictionary.keys()

                        #keys.sort(lambda a,b:cmp(int(a),int(b)))
                        #query_dictionary = extract(query_dictionary,keys[keys.index('8700'):keys.index('9000')])	

                        # if runtype == "Paironly":
                        # 	query_dictionary = extract(query_dictionary,[input_query_1,input_query_2])
                        #	pickle_data(query_dictionary,current_dir + '/pickled_data/' +"query_dictionary_" + input_query_1 + "_" + input_query_2 + ".pickle")

                else:
                        query_dictionary = load_pickle_data(current_dir + '/pickled_data/'  + pickle_filename)
               
                return  [query_dictionary,MW_test_group,isotope_masses,mods_to_check,poss_mods,mod_dict,modifications] 


##END OF getmassdict()
        

	 	
##==============================================
##FUNCTIONS
##===========================================
def readfile(name_of_file,cr = 1):
	#if screen_print:
	#	print "Loading "+name_of_file+"..."
	# print name_of_file
	# exit()
	file_data = open(name_of_file)
	file_contents = file_data.readlines()
	output = []
	for i,line in enumerate(file_contents):
                if (i <> len(file_contents)-1): # no CR on last line, apparently
                        output.append(line[:(-1* cr)]) # strip carriage returns
                else:
                        output.append(line)
	#if screen_print:
	#	print "...done"
	return output

def make_data_dictionary(data): # creates dictionary where keys are section names
	if screen_print:
		print "Creating data dictionary..."
	content = {}
	flag = 0
	for line in data:
                if "Content-Type: application" in line: #new section
                        flag = 1
                        section_name = line.split('"')[1] #key - eg. peptides, parameters
                        if "query" in section_name:
                                flag = section_name[5:] # query number
                                section_name = "query"+flag
                        continue
                elif flag > 0:
                        if flag == 1: # new section, not a query
                                content.setdefault(section_name,[]).append(line)
                        else:						
                                if line == "" or len(line) == 1:
                                        content.setdefault(section_name,[]).append("query"+str(flag))
                                else:
                                        content.setdefault(section_name,[]).append(line)
	return content


def make_query_dictionary(mass_data,peptide_data,protein_data):

	def make_peptide_dictionary(data):
		if screen_print:
			print "Creating peptide dictionary..."
		output = {}
		for line in data[1:-1]:
           #     for line in data[]:
			if "terms" not in line and "=-1" not in line and "primary" not in line and "subst" not in line and not(line is ''):
				key = int(line.split('_')[0][1:]) # query number
				pep = line.split(';')[0].split(',')
				protein_info = line.split(';')[1]
				protein_info = protein_info.replace(',','%')
				pep.append(protein_info)
				#print pep
				val = peptide_line_object(pep)
				#print val
				output.setdefault(key,[]).append(val)
		return output

	def make_MS2(data):
		queries_with_MS2_data = [] 
		
		# keys are query numbers
		# many will have the same scan numbers
		if screen_print:
			print "Creating protein dictionary..."
			
		output = {}

		if dat_file_version == "Old":
			data_length = 11
			data_range_start = 2
		elif dat_file_version == "New":
			data_length = 12
			data_range_start = 3
                for i in data:
                    if len(i) == data_length : # this is a consistency check to make sure that all the MS2 scan data is there. If not, don't bother.				
                        key = int(i[0][5:])
                        queries_with_MS2_data.append(key)
                        val = [i[1]] + i[data_range_start:data_range_start + 8]
                        output.setdefault(key,[]).append(val)
                    else:
                        print len(i)
                        print i[0]
                        print
                        print "There might be something new/different/wrong about your DAT file."
                        print
                        continue
                    #    exit()		                                
		return output, queries_with_MS2_data
		# for i in data:
		# 	if len(i) == 11 : # this is a consistency check to make sure that all the MS2 scan data is there. If not, don't bother.
                #             	data_range_start = 2
		# 	elif len(i) == 12 :
                #                 data_range_start = 3
		# 		# print len(i)
                                
		# 		# print i[0]
		# 		# print
		# 		# print "There might be something new/different/wrong about your DAT file. Aborting."
		# 		# print
		# 	#	exit()
                #         print len(i)
                #         key = int(i[0][5:])
                #         queries_with_MS2_data.append(key)
                #         val = [i[1]] + i[data_range_start:data_range_start + 8]
                #         output.setdefault(key,[]).append(val)
		# return output, queries_with_MS2_data

	output = {}
      
        #peptide_data =content["peptides"]
        #protein_data =queries
        #mass_data=content["summary"]
    peptide_dictionary = make_peptide_dictionary(peptide_data)
	peptide_queries = peptide_dictionary.keys()
	peptide_queries.sort()
	protein_dictionary, queries_with_MS2_data = make_MS2(protein_data)
	protein_queries = protein_dictionary.keys()
	for i,line in enumerate(mass_data):
		if 'qexp' in line:
			query_number = int(line.split('=')[0][4:])
			if query_number in queries_with_MS2_data:
				mass = float(line.split('=')[1].split(',')[0])
				charge = float(line.split(',')[1][0])
				qmass = mass_data[i-1]
				qintensity = mass_data[i+1]
				qmatch = mass_data[i+2]
				qplughole = mass_data[i+3]
				precursor_mass = mass * charge - 1.007825*charge
				key = str(query_number)
				val = [[query_number,mass,charge,precursor_mass,qmass,line,qmatch,qplughole,qintensity]]
				try:
					val.append(peptide_dictionary[query_number])
				except KeyError:
					val.append([])
				try:
					val.append(protein_dictionary[query_number][0])
				except KeyError:
					val.append(["No proteins"])
				#print key
				val = q_dict_object(val)
				output.setdefault(key,[]).append(val)
	return output


class q_dict_object(list):
	def __init__(self, data):
		list.__init__(self, data)
		self.charge = int(data[0][2])
		#self.precursor_MW = float(data[0][3])
		self.precursor_mz = float(data[0][1])
		self.precursor_MW = (self.precursor_mz * self.charge) - (self.charge * proton_mass)
		self.peptide_lines = data[1]
		self.peptides = [a.peptide for a in data[1]]
		self.proteins = [a.protein for a in data[1]]
		self.scores = [a.score for a in data[1]]
		if len(data[1]) > 0:
			self.top_peptide_line = data[1][0]
			self.top_peptide = data[1][0].peptide
			self.top_PME = data[1][0].PME
			self.top_mods = data[1][0].mods

			'''added tp accomodate fixed C mods - this would replace all Cs in the mod string from 0s to Cs '''
			# if 'C' in self.top_peptide:
			# 	newmods = list(self.top_mods)
			# 	for i,aa in enumerate(self.top_peptide):
			# 		if aa == 'C':
			# 			newmods[i+1] = 'C'
			# 	self.top_mods = ''.join(newmods)

			self.top_mod = data[1][0].mod
			self.top_score = data[1][0].score
			self.top_protein = data[1][0].protein
			self.top_peaks_used = data[1][0].peaks_used
			self.top_calc_mz = data[1][0].calc_mz
		else:
			self.top_peptide = "None"
			self.top_score = 0
			self.top_protein = "None"
			self.top_mods = '0'
			self.top_mod = '0'
			self.top_PME = 0
			self.top_calc_mz = 0

		self.MS2_data = data[2]
		self.MS2_ions = [[float(b) for b in a.split(':')] for a in data[2][8][6:].split(',')] #.sort(lambda x,y:cmp(x[0],y[0]))
		self.scan_number = scan_number(data)
                
class peptide_line_object(list):
	def __init__(self,data):
		list.__init__(self, data)
		self.calc_mz = float(data[1])
		self.mods = data[6]
		self.mod = int(data[6][-2])
		self.score = float(data[7])
		self.peptide = data[4]
		self.PME = float(data[2])
	#	self.score = float(data[7])
		#print self.peptide
		#print data[11]
		self.protein = [a.split(':')[0].replace('"','') for a in data[11].split('%')]
		#self.protein = [a.split('"')[1].split('|')[1] for a in data[11].split('%')]
		#self.protein = data[11].split('"')[1].split('|')[1] # yeast
		self.peaks_used = int(data[5])

def scan_number(dictionary_element):

        '''Returns scan number from a dictionary element'''
                #return int(dictionary_element[2][0].split('FinneganScanNumber')[1].split('%20')[1]) # might have to use this line depending on Mascot version
	       
        return int(dictionary_element[2][0].split('%2e')[1])





#@print_timing
def identifier(MS1,ions, isotope, cur):

	def local_isotope_check(p,i):
		#pep = 'NVTHSRMQICDQCGR' 
		#isotope = 'ARG_ARG'
		return  [sum([[mod_names[[b[0] for b in mod_names].index(a)][1] for a in i.split('_')].count(m) == peptide.count(m) for m in [a[1] for a in mod_names]])==len(mod_names) for peptide in [p]][0]

	#@print_timing #- takes about 0.01 ms
	def mods_check(peptide,com):
		# com will be list of mods... eg ([1,3])
		#print peptide,com
		aa = [modifications[m][0] for m in com] # aa = ['C','C'] or ['1','M']
		#print aa
		aa_copy = aa[:]

		'''WILL NEED TO FIX THIS FOR OTHER n-term and c-term mods. This is a one-of'''
		if '1' in aa: # n-terminal mod
			'''find the code for the n-terminal mod'''
			for mm in modifications:
				if modifications[mm][0] == '1':
					term_mods = [str(mm),'0']
					#term_mods = ['1','0']
					aa.remove('1')
		else:
			term_mods = ['0','0']

		mod_string = ['0'] * len(peptide)

		for x,p in enumerate(peptide):
			if p in aa and len(aa)>0:
				mod_string[x] = str(com[aa_copy.index(p)])
				aa.remove(p)
		if len(aa) > 0:
			return False
		else:
			mod_string = ''.join(mod_string)
			return term_mods[0] + mod_string + term_mods[1]


	def getSubset(l,h):
		min_idx=bisect.bisect_left(tryptic_masses, l)
		max_idx=bisect.bisect_right(tryptic_masses, h, min_idx)
		peptides_to_check = [[mass,tryptic_dictionary[mass]] for mass in tryptic_masses[min_idx:max_idx]]
		output = []
		for line in peptides_to_check:
			#print line
			mass = line[0]
			peptides = line[1]
			peptides.sort(lambda a,b:cmp(a[0],b[0]))
			p = []
			for peptide in peptides:
				pep = peptide[0]
				pro = peptide[1]
				if pep in p:
					output[-1][2]+="|"+pro
				else:
					p.append(pep)
					output.append([mass,pep,pro])

		peptides_to_check = output
		return peptides_to_check

	#@print_timing
	def get_peptides_to_test(C,MS1):
		peps_to_test = []
		for com in C:
			mass = sum([modifications[a][1] for a in com])
			m = MS1 - mass
			Da = (identifier_width_ppm * m) / 1000000.0
			l = m - Da
			h = m + Da
			#print "Getting around mass",m
			peps = getSubset(l,h)
			#for line in peps:
			#	print line
			# print
			# ps = [a[1] for a in peps]
			# ps.sort()
			peps = [a for a in peps if local_isotope_check(a[1],isotope)] # removes peptides that are not consistent with the isotope
			#new_peps = []
			for pep in peps:
				#print pep # [1874.86934, 'EAESSGLPWSCQVGGRGR', 'IPI00791344']
				mods = mods_check(pep[1],com) # will check the peptide for compatibilty with the mods. If OK, will return a possible mod string

				if mods:
					#new_peps.append([a for a in pep] + [mods])
					peps_to_test.append([a for a in pep] + [mods] + [com])
		return peps_to_test


	if screen_print:
		print 
		print "*" * 80
		print
		print "Start of IDENTIFIER"
		print
		print "*" * 80
		print

	''' Attempt at adding fixed mods - while this would significantly cut down on search time, it would be difficult to test every possibility'''
	# if 'fixed' not in modifications:
	# 	'''Assuming: {0: ['N', 0.0], 1: ['C', 57.021469], 2: ['1', 27.994919], 3: ['M', 15.994919], 4: ['K', 8.02684], 5: ['R', 10.00827]}'''
	# 	combos = ['1','11','3','33','13','113','133','1133','2','21','211','23','233','213','2113','2133','21133']
	# else:
	# 	'''Assuming: {0: ['N', 0.0], 1: ['1', 27.994915], 2: ['M', 15.994915], 3: ['K', 8.014199], 4: ['R', 10.008269], 'fixed': [['C', 57.021464]]}'''
	# 	combos = ['2','22','12','122']


	combos = ['1','11','3','33','13','113','133','1133','2','21','211','23','233','213','2113','2133','21133']

	output = []
	for c in combos:
		output.append([int(a) for a in c])
	combos = output
	combos = [[''],combos]


	all_peps = []
	all_random_peps = []
	all_PMEs = []

	for Z,C in enumerate(combos):
		peps = get_peptides_to_test(C,MS1)
		PMEs = [calc_identifier_PME(line[1] + '_' + line[3],MS1) for line in peps]
		all_PMEs += PMEs
		all_peps += peps
		#all_random_peps += random_peps


	if len(all_peps) > 0:
		all_random_peps = []
		for i in range(bootstrap_n):
			pep = choice(all_peps)
			p = list([a for a in pep[1][:-1]])
			shuffle(p)
			random_pep = p + [pep[1][-1]]
			random_pep = ''.join(random_pep)
			com = pep[4]
			random_mods = mods_check(random_pep,com)
			if random_mods:
				all_random_peps.append([random_pep] + [random_mods])



	if len(all_peps) == 0:
		if screen_print:
			print "Found no peptides in that range"
		winning_peptide = 'None'
		winning_protein = 'None'
		scores = [0]
		peps_and_scores = []
		b = 0
		B = 0
		V = 0
		t = 0
		t_bar = 0
	else:
		if screen_print:
			print
			print "About to test these peptides around mass:",MS1
			for pep in all_peps:
				print pep[1],pep[3]

		output = [[],[]]
		# print "LENGTH",len(all_peps),len(all_random_peps)
		random_flag = [False,True]
		for i,P in enumerate([all_peps,all_random_peps]):
			for line in P:
				if i == 0:
					peptide = line[1]
					#proteins = line[2]
					mods = line[3]
					if screen_print:
						print "Now testing",peptide
						print					
				else:
					peptide = line[0]
					mods = line[1]

				mod_string = create_mod_string(mods[1:-1],[mods[0],mods[-1]], peptide)
				fragment_ions = fragment(peptide,mod_string,random_flag[i]) # random flag used to suppress printing of fragmentation table for random peptides.					
				ions_,total_consecutive,[scoreb,scorey],percent_intensity = score_peptide(fragment_ions, ions, peptide,random_flag[i])

				if screen_print and i==0: # don't do this for random peptides
					#print "Scores for peptide "+peptide+" = ",scoreb,scorey,scoreb+2.0*scorey, "[score b, score y, score b + 2 x score y]"
					print "Scores for peptide "+peptide+" = ",scoreb,scorey,total_consecutive,scoreb+2.0*scorey + 5.0 * float(total_consecutive), "[score b, score y, score b + 2 x score y + 5*totalconsecutive]"
					print
					print "############"
					print
					print
				output[i].append([peptide, total_consecutive,[scoreb,scorey],percent_intensity])
				# print "b",scoreb
				# print "y", scorey
				# print "consec",total_consecutive
				# print "score",scoreb+2.0*scorey + 5.0 * float(total_consecutive)
				# exit()
		#scores_peps = [a[2][0] + 2.0 * a[2][1] for a in output[0]] # simple score of B +  2Y
		scores_peps = [a[2][0] + 2.0 * a[2][1] + 5.0 * float(a[1]) for a in output[0]] # simple score of B +  2Y

		scores_random_peps = [a[2][0] + 2.0 * a[2][1] for a in output[1]]
		#scores = zip(scores,PMEs, [(float(a[0]) / (0.5 + abs(float(a[1])) ) ) for a in zip(scores, PMEs)]        )
		scores = zip(scores_peps,all_PMEs)
		peps_and_scores = zip(all_peps,scores)


		t = average(scores_peps) # this is 't'

		t_bar = average(scores_random_peps) # this is t-bar
		#st = stdev(scores_random_peps)
		R = len(scores_random_peps)
		B = 1.0/float(R) * sum([t_star-t for t_star in scores_random_peps])
		V = 1.0/(float(R)-1.0) * sum([((t_star - t_bar) * (t_star - t_bar)) for t_star in scores_random_peps])

		# print "t",t
		# print "t_bar",t_bar
		# #print "stdev",st
		# print "R",R
		# print "B", B
		# print "V",V

		peps_and_scores.sort(lambda a,b:cmp(abs(a[1][1]),abs(b[1][1])))# sort by ascending absolute PME
		# for line in peps_and_scores:
		# 	print line
		# print
		# print
		if FDR == 0.05:		
			zstat = 1.96 #95%
		elif FDR == 0.01:
			zstat = 2.58 # 99%
		elif FDR == 0.001:
			zstat = 3.29 # 99.9%
		elif FDR == 0.0001:
			zstat = 3.89 # 99.9%
		else:
			print "FDR Undefined - using 5%"
			zstat = 1.96


		#a = (t - B) - (math.sqrt(V) * zstat)
		b = (t - B) + (math.sqrt(V) * zstat)

		# if len(peps_and_scores) > 0 and len(scores_random_peps)>0 and V>0.0:
		# 
		# 	for line in peps_and_scores:
		# 		peptide = line[0][1]
		# 		score = line[1][0]
		# 		PME = line[1][1]
		# 		z = stats.lz(scores_random_peps,score)
		# 		l = stats.lzprob(z)


			#pyplot.hist(scores_random_peps)
			# counter = 0
			# for line in peps_and_scores:
			# 	pyplot.axvline(x=line[1][0], color='r')
				# if counter < 10:
				# 	print line[0][1],line[1][0],line[1][1]
				# 	counter += 1
			# pyplot.figtext(0.5,0.8,str(avg))
			# pyplot.figtext(0.5,0.6,str(std))
			# pyplot.title(winning_peptide)
			#pyplot.show()
			#pyplot.clf()
	# for line in peps_and_scores:
	# 	print line
	# print
	# print		
	master_list = peps_and_scores

	if screen_print3:
		print "CUTO FF SCORE",b
		O = []
		if len(master_list) > 0:
			for line in master_list:
				winning_peptide = line[0][1]
				winning_proteins = line[0][2]
				winning_mods = line[0][3]
				winning_peptide = winning_peptide + '_' + winning_mods
				score = line[1][0]
				PME = line[1][1]
				newline = [winning_peptide, score, PME]
				O.append(newline)
			out = sys.stdout
			pprint_table(out, O)
			print
		else:
			print "NO OUTPUT"



	#master_list.sort(lambda a,b:cmp(b[1][0],a[1][0])) # sort by scores
	winning_peptide = "None"
	winning_proteins = "None"
	high_score = 0
	winner_PME = 0

	# first try unmodified
	unmodified = [line for line in master_list if int(line[0][3]) == 0]
	modified = [line for line in master_list if int(line[0][3]) > 0]

	if screen_print3:
		for line in unmodified:
			print line
		print
		for line in modified:
			print line

		print
		print '___'
	'''This statement sorts by PME but groups all similar (but not exact PMEs together). So, when dividing by 1000, it groups PMEs the same in the 1000th place.'''
	unmodified.sort(lambda a,b:(cmp(abs(math.trunc(a[1][1]*PME_rounder)/PME_rounder),abs(math.trunc(b[1][1]*PME_rounder)/PME_rounder)) or cmp(b[1][0],a[1][0])))
	modified.sort(lambda a,b:(cmp(abs(math.trunc(a[1][1]*PME_rounder)/PME_rounder),abs(math.trunc(b[1][1]*PME_rounder)/PME_rounder)) or cmp(b[1][0],a[1][0])))

	master_list.sort(lambda a,b:(cmp(abs(math.trunc(a[1][1]*PME_rounder)/PME_rounder),abs(math.trunc(b[1][1]*PME_rounder)/PME_rounder)) or cmp(b[1][0],a[1][0])))


	if screen_print3:
		for line in unmodified:
			print line
		print
		print '___'
		for line in modified:
			print line
			print line
		print
		print '___'




	# for line in unmodified:
	# 	print line
	# print
	# print
	# for line in modified:
	# 	print line
	# print
	# exit()

	flag = 0
	for PME_cutoff in PME_cutoffs: # this should search in this order. Unmodified 0-10ppm
		#for samples in [unmodified, modified]:
		for samples in [master_list]:	
			for line in samples:
				score = line[1][0]
				PME = line[1][1]
				# if abs(PME)<10:
				# 	print line
				#print line
				if score > b and (abs(PME) >= PME_cutoff[0] and abs(PME) <=PME_cutoff[1]):
					# winner
					winning_peptide = line[0][1]
					winning_proteins = line[0][2]
					winning_mods = line[0][3]
					winning_peptide = winning_peptide + '_' + winning_mods
					high_score = score
					winner_PME = PME
					flag = 1
					break
			if flag == 1:
				break
		if flag == 1:
			break

	bin_size = len(master_list)
	R = 0
	runner_up = 'None'

	# if len(master_list)>0:
	# 	winning_peptide = master_list[0][0][1]
	# 	winning_proteins = master_list[0][0][2]
	# 	winning_mods = master_list[0][0][3]
	# 	winning_peptide = winning_peptide + '_' + winning_mods
	# 	high_score = master_list[0][1][0]
	# 
	# 	if len(master_list) >1:
	# 		runner_up = master_list[1][0][1] + "_" + str(master_list[1][1][0])
	# 		bin_size = len(master_list) 
	# 		
	# 		if master_list[0][1][0] > 0 and master_list[1][1][0] > 0:
	# 			R = master_list[0][1][0]/master_list[1][1][0] # the default.
	# 		elif master_list[0][1][0] > 0:
	# 			R = "" # first one is >0 but rest are zero
	# 		else:
	# 			R = 0 # all are zero
	# 	else: # len = 1
	# 		R = ""
	# 		runner_up = 'None'
	# 		bin_size = len(master_list)
	# 		
	# else:
	# 	R = 0
	# 	winning_peptide = 'None'
	# 	winning_proteins = 'None'
	# 	high_score = 0
	# 	runner_up = 'None'
	# 	bin_size = len(master_list)
	# 
	# 
	# print 'Winning score',high_score, winning_peptide
	# print
	# 
	# if screen_print:
	# 	print
	# 	print "After all testing, R=",R
	# 	print
	# 	O = []
	# 	for line in master_list:
	# 		#print line
	# 		O.append([line[0][1],line[0][3],line[1]])
	# 	out = sys.stdout
	# 	pprint_table(out, O)
	# 	print
	# 
	# if screen_print:
	# 	print
	# 	print
	# 	print "Overall winning peptide, winning_protein, cutoff, score = ",winning_peptide, winning_proteins, R, high_score, runner_up, bin_size
	# 	print
	# 	print 
	# 	print "*" * 80
	# 	print
	# 	print "End of IDENTIFIER"
	# 	print
	# 	print "*" * 80
	# 	print
	if screen_print3:
		print winning_peptide, winner_PME, high_score,b
		print
	return winning_peptide, winning_proteins, b, high_score, runner_up, bin_size, [B,V,t,t_bar]

def identify_peptides(matches_3,Q):
	# each line is ([[q1,q2],c,[iso,isotopologue,PME], final_score])
	if use_postgres:
		#conn = psycopg2.connect("host=proteomics3.ath.cx user=postgres password=postgrespasswd dbname=tryptic")
		conn = psycopg2.connect("user=volcs0 dbname=tryptic")
		cur = conn.cursor()
	else:
		cur = None

	#@print_timing
	def create_tryptic_dictionary():
		#data = readfile(current_dir + '/databases/' + "tryptic_IPI_postgres.txt")
		data = readfile('/Volumes/SpeedDisk/tryptic_IPI_postgres.txt')
		#print len(data)
		output = {}
		for line in data:
			line = line.split('\t')
			mass = float(line[0][:-9])
			peptide = line[1]
			protein = line[2]#[:-1]
			#print mass, peptide, protein

			output.setdefault(mass,[]).append([peptide,protein])
		return output

	if not use_postgres:	
		print "Loading tryptic dictionary"
		global tryptic_masses, tryptic_dictionary
		tryptic_dictionary = create_tryptic_dictionary()
		tryptic_masses = tryptic_dictionary.keys()
		tryptic_masses.sort()

	# global tryptic_dictionary, tryptic_masses
	# tryptic_dictionary = create_tryptic_dictionary()
	# tryptic_masses = tryptic_dictionary.keys()
	# tryptic_masses.sort()	

	# global tryptic_masses, tryptic_dictionary
	# tryptic_dictionary = load_pickle_data(current_dir + "/" + "tryptic_IPI_small2.pickled")
	# tryptic_masses = tryptic_dictionary.keys()
	# tryptic_masses.sort()


	#print "Identifying peptides"


	for X,line in enumerate(matches_3):
		isotopologue = line[2][1]

		if not real_time:
			if X%100 == 0: print X, len(matches_3)

		q1 = line[0][0]
		q2 = line[0][1]
		#print "QUERIES",q1,q2
		# print "*" * 80
		# print "QUERY",q1,q2

		D1 = Q[q1][0]
		D2 = Q[q2][0]
		# print
		# print "MASCOT",D1.top_peptide
		# print "MASCOT",D2.top_peptide
		# print
		if screen_print3:
			print "-----"
			if D1.top_peptide == D2.top_peptide:
				print "MASCOT match", D1.top_peptide, q1, q2

		iso = line[2][0]
		ionsL = [a[0] for a in line[1][0]],[a[0] for a in line[1][1]],[a[0] for a in line[1][2]],line[1][4][0],line[1][5][0] # non-shifers, shifters, half-shifters, sum_intensities, max intensity
		t1 = ionsL[0]
		t1 = unique(t1)
		t1.sort()
		t2 = ionsL[1]
		t2 = unique(t2)
		t2.sort()
		t3 = ionsL[2]
		t3 = unique(t3)
		t3.sort
		ionsL = (t1,t2,t3,ionsL[3],ionsL[4])
		ionsH = [a[1] for a in line[1][0]],[a[1] for a in line[1][1]],[a[1] for a in line[1][2]],line[1][4][1],line[1][5][1]
		t1 = ionsH[0]
		t1 = unique(t1)
		t1.sort()
		t2 = ionsH[1]
		t2 = unique(t2)
		t2.sort()
		t3 = ionsH[2]
		t3 = unique(t3)
		t3.sort()
		ionsH = (t1,t2,t3,ionsH[3],ionsH[4])


		ions = [ionsL, ionsH]



		#print "PEPTIDE",query_dictionary[q1][0].top_peptide
		precursor_mw_L = query_dictionary[q1][0].precursor_MW

		''' ADDED 3/10/11'''
		if isotopologue == "L":
			precursor_mw_L -= first_isotopologue
		elif isotopologue == "2L":
			precursor_mw_L -= second_isotopologue


		#precursor_mw_H = query_dictionary[q2][0].precursor_MW

		#peptideL = identifier(precursor_mw_L,ionsL, iso, cur)
		#peptideH = identifier(precursor_mw_H,ionsH, cur)


		winning_peptide, winning_protein, cutoff_score, high_score, runner_up, bin_size, [bias, variance, av_peptides, av_random_peptides] = identifier(precursor_mw_L,ionsL, iso, cur)
		line.append([winning_peptide, winning_protein, cutoff_score, high_score, runner_up, bin_size, [bias, variance, av_peptides, av_random_peptides]])

	return matches_3



def extract(d,keys):

	# used to create a smaller dictionary subset for testing
	return dict((k, d[k]) for k in keys if k in d)

def pickle_data(data, filename):
	newfile = open(filename,"w")
	cPickle.dump(data,newfile)
	newfile.close()
	return
###test code#####
