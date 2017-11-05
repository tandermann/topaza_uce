import os
import re
import csv
import argparse


# Complete path function
class CompletePath(argparse.Action):
	"""give the full path of an input file/folder"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def get_args():
	parser = argparse.ArgumentParser(
		description="Get the distance to center of each extracted SNP",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		"--logfile",
		required=True,
		action=CompletePath,
		help="The screen output from the SNP extraction script as a text file"
	)
	parser.add_argument(
		"--length_info",
		required=True,
		action=CompletePath,
		help="supply a file that contains the locus name and the length of the locus (each in a separate line) for all alignments of interest"
	)
	parser.add_argument(
		"--output",
		required=True,
		action=CompletePath,
		help="where to store the results (dir)"
	)
	parser.add_argument(
		'--multi_snp',
		action='store_true',
		default=False,
		help='Use this flag if you extracted all snps (not only one per locus)'
	)
	return parser.parse_args()

args = get_args()
snp_logfile = args.logfile
length_log = args.length_info
outdir = args.output


# get extracted positions

logfile = open(snp_logfile,'r')
alignments = []
extracted_positions = []
if args.multi_snp:
	for line in logfile:
		if line.startswith('processing file:'):
			alignment = line.split(':')[1]
			alignment =  re.sub('\s+', '', alignment)
			alignments.append(alignment)
		elif line.startswith("['"):
			pass
		elif line.startswith("[]"):
			pass
		elif line.startswith("["): 
			line = line.rstrip()
			pos = re.sub('[\[\]]','',line)
			extracted_positions.append(pos)
		elif line.startswith('no SNP extraction performed'):
			pos = 'NA'
			extracted_positions.append(pos)
else:	
	for line in logfile:
		if line.startswith('processing file:'):
			alignment = line.split(':')[1]
			alignment =  re.sub('\s+', '', alignment)
			alignments.append(alignment)
		elif line.startswith('sampling position'):
			pos = line.split(' ')[-1]
			pos = pos.rstrip()
			extracted_positions.append(pos)
		elif line.startswith('no SNP extraction performed'):
			pos = 'NA'
			extracted_positions.append(pos)
#make a dicitonary from the two lists
al_pos_dict = {}
for al in alignments:
	index = alignments.index(al)
	al_pos_dict.setdefault(al,extracted_positions[index])



# get alignment lengths

lengthfile = open(length_log,'r')
l_alignments = []
length_list = []
for line in lengthfile:
	if line.startswith('uce-'):
		l_alignments.append(line.rstrip())
	else:
		length_list.append(line.rstrip())
#make a dicitonary form the two lists
al_length_dict = {}
for al in l_alignments:
	index = l_alignments.index(al)
	al_length_dict.setdefault(al,length_list[index])



# now calculate the distance to the center for each variable position
dist_count_dict = {}
for al in al_length_dict:
	length = al_length_dict[al]
	extracted = al_pos_dict[al]
	center = int(length)/2
	var_pos_list = al_pos_dict[al].split(', ')
	for pos in var_pos_list:
		dist = ''
		if not pos == 'NA':
			pos = int(pos)
			if pos > center:
				dist = pos-center
			elif pos < center:
				dist = center-pos
			else:
				dist = 0
		else:
			continue
		if dist in dist_count_dict:
			dist_count_dict[dist] += 1
		else:
			dist_count_dict.setdefault(dist,1)

dist_count = open("%s/distance_to_center_count.txt" %outdir, "wb")
dist_count_log=csv.writer(dist_count, delimiter='\t')
for key in dist_count_dict:
	dist_count_log.writerow([key,dist_count_dict[key]])
#input file = log file from snp extraction script
#alignment-varpos-dict = transform file into dict, with name of locus as key and extracted variable positions as value(s)
#alignment-length-dict = transform file into dict, every first line is key and every second is value
#distance_to_center-count-dict = {}
#for alignment, length in alignment-lenth-dict:
#	center = length/2
#	varpos = alignment-varpos-dict[alignment]
# 	for position in varpos:
#		dist = ''
#		if position > center:
#			dist = position - center
# 		elif position < center: 
#			dist = center-position
# 		else:
#			dist = 0
#		if dist in distance_to_center-count-dict:
#			distance_to_center-count-dict[dist] +1
#		else:
#			distance_to_center-count-dict.setfault(dist,1)
# write distance_to_center-count-dict into tab delimited file
# key /t value