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
from padnums import pprint_table
##other py files
import mascotparser

#from functionfile import *



#import stats
from PeptideFragmentSingleton import PeptideFragment
## from Gene Selkov
# '''uncomment this line to allow HTML reporting of errors when running as a webservice'''
#cgitb.enable()

# global f
# f = PeptideFragment(1)

#exit()



##======================
##Main Funciton
##======================


# global f
# f = PeptideFragment(1)
global current_dir
#current_dir="/home/lteng/Desktop/Projects/6MassSpec/ling_code"
global dat_file_version
dat_file_version = "Old"

global screen_print, screen_print2, screen_print3
global use_postgres
global single_run # used for batching 
single_run = True
global real_time # experimental - used for real-time peptide ID (eg. finds a pair, identified it, and so on)
real_time = False
#single_run_file_folder = "data"
single_run_file_folder = ""
form = cgi.FieldStorage() 
global single_run_filename, pickled_data, pickle_filename
single_run_filename = form.getvalue('datafile')
if single_run_filename == None:
    # this is dat_file_version = 'New'
    #single_run_filename = "F001524.dat"
    single_run_filename = "/Volumes/SpeedDisk/1561.dat"
    # this is dat_file_version = 'Old'
    #single_run_filename = "Mann_05_K8R10_F2901.dat"


pickle_filename = form.getvalue('pickled_datafile')
if pickle_filename == None:
    pickle_filename = "None"
global input_query_1, input_query_2
if pickle_filename <> "None" and pickle_filename <> None:
    pickled_data = True
    if "small" in pickle_filename:
        input_query_1 = input_query_2 = 0
    else:
        input_query_1 = pickle_filename.split('_')[2]
        input_query_2 = pickle_filename.split('_')[3].split('.')[0]
else:
    pickled_data = False


global runtype
runtype = form.getvalue('runtype')
if runtype == None:
    runtype = "FullRun"
# 
runtype = "Paironly"
input_query_1 = "3106"
input_query_2 = "3413"
use_postgres = False

if runtype == 'FullRun':	
    print "Content-type: text/html"
    print
    print "<html><head>"
    print ""
    print "</head><body></body>"
    use_postgres = False
    input_query_1 = "N/A"
    input_query_2 = "N/A"

# elif runtype == "Paironly":
# 	print "Content-Type: text/plain;charset=utf-8"
# 	print
# 	use_postgres = True
# 	if not pickled_data:
# 		input_query_1 = form.getvalue('query_1')
# 		input_query_2  = form.getvalue('query_2')

output = form.getvalue('output')


if output == "screen":
    screen_print2 = True
    screen_print = True
    write_output = False
elif output == "file" or output == None:
    screen_print2 = False
    screen_print = False
    write_output = True
elif output == "both":
    screen_print2 = True
    screen_print = True
    write_output = True
else:
    print output
    exit()

print_flag = True	
screen_print = print_flag
screen_print3 = print_flag	
screen_print2 = print_flag

global peak_cutoff
peak_cutoff = form.getvalue('peak_cutoff')
if peak_cutoff == None:
    peak_cutoff = 0.05

global MS1_cutoff
MS1_cutoff = form.getvalue('precursor_mass_error')
if MS1_cutoff == None:
    MS1_cutoff = 10.

global MS2_cutoff
MS2_cutoff = form.getvalue('fragmentation_mass_error')
if MS2_cutoff == None:
    MS2_cutoff = 1500.

global identifier_width_ppm
identifier_width_ppm = form.getvalue('identifier_width_ppm')
if identifier_width_ppm == None:
    identifier_width_ppm = 50.


global max_peak_difference
max_peak_difference = form.getvalue('max_peak_difference')
if max_peak_difference == None:
    max_peak_difference = 0.50
# print "<br>"
# print "############################"
# print "<br>"
# print "Here are the run parameters.<br>"
# print "Filename:",single_run_filename,"<br>"
# print "Pickle filename:",pickle_filename,"<br>" 
# print "Query 1",input_query_1,"<br>"
# print "Query 2",input_query_2,"<br>"
# print "MS1 cutoff (ppm)",MS1_cutoff,"<br>"
# print "MS2 cutoff (ppm)",MS2_cutoff,"<br>"
# print "Peak cutoff",peak_cutoff,"<br>"
# print "Max peak difference",max_peak_difference,"<br>"
# print "<br>"
# print "################################<br>"
# print "<br>"

global proton_mass
global fignum
global FDR
global PME_cutoffs
global PME_rounder
global bootstrap_n
bootstrap_n = 1000 
#PME_cutoffs = [[0.0,1.0],[1.0,2.0],[2.0,3.0],[3.0,10.0],[10.0,50.0]]
#  used to direct the peptide matching. So, identifier will iterate through the first PME range
#  (both mod and unmod) before going to the next range
PME_cutoffs = [[0.0,10.0],[10.0,50.0]] #
FDR = 0.01 # z-value cutoff. So, 0.05 would be 5% FDR
PME_rounder = 0.2 # when sorting peptide possibilities for Identifier, decides how to bunch them together. e.g. 1000 means sort all same at the 1000th place. Use 1 to sort at the ones place. 0.2 will sort every 5, and so on
fignum = 1
proton_mass = 1.00727646688
neutron_mass = 1.0086649156
first_isotopologue = 1.0033548378 # http://en.wikipedia.org/wiki/Isotopes_of_carbon
second_isotopologue  = 2 * first_isotopologue
ppm_cutoff = int(MS1_cutoff) # ppm
scan_width = 100
ms2_ppm_tolerance = int(MS2_cutoff)
min_ion_percent_intensity = float(peak_cutoff)# mininum to keep
max_ion_difference = float(max_peak_difference) # max allowable difference between ion percentages to keep



#print "FINISHED"

# ##======================
# ## Main Funciton
# ##======================
#main():
#global query_dictionary
global modifications
global mods_to_check
global poss_mods
global MW_test_group
global mod_names
global isotope_masses
global mod_dict
global current_dir
global runtime
now = datetime.datetime.now()
runtime = now.strftime("%Y%m%d%H%M%S")


mod_names = [["LYS",'K'],["ARG","R"]]
   


filename="/home/lteng/Desktop/Projects/6MassSpec/ling_code/databases/Mascot_file_Mann_HelaEpo_05_F2893.dat"
parse_result = mascotparser.getmassdict(filename)
Q =parse_result[0] #[[query_number,mass,charge,precursor_mass,qmass,line,qmatch,qplughole,qintensity]]
MW_test_group =parse_result[1]
isotope_masses = parse_result[2]
mods_to_check = parse_result[3]
poss_mods = parse_result[4]
mod_dict = parse_result[5]
modifications = parse_result[6]
# pairs_1 = validator_1(Q,mods_to_check,poss_mods, mod_dict,modifications) 
# queries_in_pairs_1 = unique(flatten([[a.split(',')[0],a.split(',')[1]] for a in pairs_1[1:]]))
# queries_in_pairs_1.sort()
# testresult =unpaired_analysis(queries_in_pairs_1,"Unpaired_analysis_Validator_1.xls")



def test():
    validator_1(Q,mods_to_check,poss_mods, mod_dict,modifications) 
  
    return 
    
    
#valid01res = validator_1(query_dictionary)
#Q=Q
################
#FUNCTIONS
################
def validator_1(Q,mods_to_check,poss_mods, mod_dict,modifications):
    
    if screen_print:
        print
        print "##################################"
        print "#########   VALIDATOR 1 ##########"
        print "##################################"
        print
    queries = Q.keys() #[query_number,mass,charge,precursor_mass,qmass,line,qmatch,qplughole,qintensity]
    queries.sort(lambda x,y: cmp(Q[x][0].scan_number,Q[y][0].scan_number))
    scan_numbers = [Q[a][0].scan_number for a in queries]
    D = zip(scan_numbers,queries)
    pairs = []

    for i,line in enumerate(D):
        q1 = line[1]
        c1 = Q[q1][0].charge
        pep1 = Q[q1][0].top_peptide
        s1 = line[0]
        ###
        ### Find the scan window 
        ###
        low = max(s1-scan_width,0)
        while low not in scan_numbers and low <=s1:
            low += 1
        high = min(s1+scan_width, scan_numbers[-1])
        while high not in scan_numbers and high >= s1:
            high -= 1
        i1 = scan_numbers.index(low)
        i2 = scan_numbers.index(high)
        R2 = D[i1:i2]
        #
        #
        #print mods_to_check		
        for line2 in R2: # iterate through the scans
            q2 = line2[1]
            mz1 = Q[q1][0].precursor_mz
            mz2 = Q[q2][0].precursor_mz
            pep2 = Q[q2][0].top_peptide
            c2 = Q[q2][0].charge
            s2 = line2[0]

            if q1 <> q2 and pep1 == pep2 and c1 == c2 and mz1<>mz2:
                if mz1<mz2:
                    pepL = pep1
                    pepH = pep2
                    modsL = Q[q1][0].top_mods[1:-1]
                    modsH = Q[q2][0].top_mods[1:-1]
                    qL = q1
                    qH = q2
                else:
                    pepL = pep2
                    pepH = pep1
                    modsL = Q[q2][0].top_mods[1:-1]
                    modsH = Q[q1][0].top_mods[1:-1]
                    qL = q2
                    qH = q1

                if modsL <> modsH and len(pepL) == len(modsL) and len(pepH) == len(modsH) and mods_consistent(mods_to_check,mod_dict,pepL,modsL,'L') and mods_consistent(mods_to_check,mod_dict,pepH,modsH,'H') and [qL,qH] not in pairs: 
                    pairs.append([qL,qH])
#Q:[query_number,mass,charge,precursor_mass,qmass,line,qmatch,qplughole,qintensity]
    pairs.sort(lambda a,b:cmp(int(Q[a[0]][0].scan_number),int(Q[b[0]][0].scan_number)))
    output = []
    line1 = "Query 1, Query 2, Peptide, Mods 1, Mods 2, Scan 1, Scan 2, Score 1, Score 2, PME, ? Isotopologue, Proteins"
    output.append(line1)
    for line in pairs:
       #line =pairs[0]	
        q1 = line[0]
        q2 = line[1]		
        mods1 = Q[q1][0].top_mods[1:-1]
        mods2 = Q[q2][0].top_mods[1:-1]
        dif_modsL = [a[0] for a in zip(mods1,mods2) if a[0] <> a[1]] # find the differences between the two
        print "0"
        dif_modsH = [a[1] for a in zip(mods1,mods2) if a[0] <> a[1]] # find the differences between the two 
        subL = sum([modifications[int(mod)][1] for mod in dif_modsL])
        subH = sum([modifications[int(mod)][1] for mod in dif_modsH])
        sub = subH-subL
        ME = (Q[q2][0].precursor_MW - sub) - Q[q1][0].precursor_MW 
        PME = ME - neutron_mass * (ME > 0.5) + neutron_mass * ( ME < -0.5 )
        print "1"
        if ME > 0.5 or ME < -0.5:
            add = "**"
        else:
            add = ""
        print "2"
        ppm_calc = PME/Q[q1][0].precursor_MW *1000000
        proteins = Q[q1][0].top_protein
        prots = ",".join(proteins)		
        score1 = Q[q1][0].top_score
        score2 = Q[q2][0].top_score		
        newline =  [q1,q2,Q[q1][0].top_peptide,"'"+str(mods1),"'"+str(mods2),Q[q1][0].scan_number,Q[q2][0].scan_number,score1, score2,abs(ppm_calc),add, prots]
        newline = ",".join([str(a) for a in newline])
        output.append(newline)
    print "Writing output to csv file "
    writefile(output,"Val_1_pairs_" + runtime + ".csv")
    queries_in_pairs_1 = unique(flatten([[a.split(',')[0],a.split(',')[1]] for a in output[1:]]))
    queries_in_pairs_1.sort()
    print "Doing unpaired_analysis..."
    unpaired_analysis(Q,queries_in_pairs_1,"Unpaired_analysis_Validator_1.xls")
    print "Validator_1 is Done."
    return 

##
def mods_consistent(mods_to_check, mod_dict,pep,mods,kind):

    M = ['K','R']

    for l in zip(pep,mods):
        p = l[0]
        m = int(l[1])
        if p in M and kind == 'L' and m <> 0:
            return False
        elif p in M and kind == 'H' and m <> mod_dict[p]:
            return False
        elif p not in M and m in mods_to_check:
            return False
    return True
	
#
def writefile(output,filename):
    newfile = open(filename, "w")
    for a in output:
        newfile.write(a+"\n")
    newfile.close()
    return
##
def unique(s):
    n = len(s)
    if n == 0:
        return []
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u
##
def unpaired_analysis(Q,queries_in_pairs, filename):
    # [['7733', '7753'], ['ARG', 10.00827], 'N', -0.12918151875875858]
    #queries_in_pairs = unique(flatten([ [int(b) for b in a[0]] for a in pairs]))
    #queries_in_pairs.sort()

    queries = Q.keys()

    unpaired = [int(query) for query in queries if int(query) not in queries_in_pairs]
    unpaired = unique(unpaired)
    unpaired.sort()
    header = ['Query','Scan','mz','charge','mw','Peptide','Mods','Score','Protein']
    output = [header]
    for q in unpaired:
        D = Q[str(q)][0]
        scan = D.scan_number
        mz = D.precursor_mz
        charge = D.charge
        mw = D.precursor_MW
        peptide = D.top_peptide
        mods = D.top_mods
        score = D.top_score
        protein = '|'.join(D.top_protein)
        newline = [q,scan,mz,charge,mw,peptide,mods,score,protein]
        output.append(newline)

    write_excel_output(output,filename)
	
##
def flatten(x):
    """flatten(sequence) -> list"""
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result
##
def write_excel_output(data,filename):
    '''generic routine to write an excel table'''


    #header = data[0]
    #data = data[1:]

    wb = Workbook()
    ws_unpaired_data = wb.add_sheet('unpaired_data',cell_overwrite_ok=True)
    ws_unpaired_data.panes_frozen = True
    ws_unpaired_data.remove_splits = True
    ws_unpaired_data.horz_split_pos = 1
    ws_unpaired_data.horz_split_first_visible = 2

    style_none = easyxf('font: name Courier;')
    style = easyxf('font: name Courier; align: horizontal center;')
    dec_style = easyxf('font: name Courier;', num_format_str='0.000')

    for x,line in enumerate(data):
        for y,l in enumerate(line):
            #print l	
            ws_unpaired_data.row(x).write(y,l, style)	

    wb.save(filename)
