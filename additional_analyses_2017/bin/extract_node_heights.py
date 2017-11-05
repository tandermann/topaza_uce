import os
import re
import csv
#import dendropy
from operator import itemgetter

#treefile = raw_input("Give path to tree file: ")
treefile = "/Users/tobias/Desktop/figures_topaza_review/stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/species.trees"
burnin = 0.1
process_trees = 'all'
clade_list = []
if 'simulated_contigs' in treefile:
	clade_list = ['((d1,d2),(e1,e2))','(y1,(z1,z2))','((x1,x2),(y1,(z1,z2)))', '(((d1,d2),(e1,e2)),((x1,x2),(y1,(z1,z2))))']
elif 'simulated_alleles' in treefile:
	clade_list = ['((e1,e2),(d1,d2))','((z1,z2),y1)','(((z1,z2),y1),(x1,x2))', '(((e1,e2),(d1,d2)),(((z1,z2),y1),(x1,x2)))']
elif 'simulated_allele_snps' in treefile:
	clade_list = ['((E1_,E2_),(D1_,D2_))','((Z1_,Z2_),Y1_)','(((Z1_,Z2_),Y1_),(X1_,X2_))','((((Z1_,Z2_),Y1_),(X1_,X2_)),((E1_,E2_),(D1_,D2_)))']
#######################################################################################################


def get_metadata(treefile):
	reader = csv.reader(open(treefile, 'r'), delimiter=',')
	reader = list(reader)
	header_index = ''
	for row in reader:
		if len(row) > 0 and row[0].startswith('tree STATE_0'):
			header_index = int(reader.index(row)-1)
			print 'skipping first', header_index, 'lines (header)'
		
	info = reader[0:header_index]
	for row in info:
		if len(row) > 0 and row[0].strip() == 'Translate':
			begin_translate = info.index(row)
		
	translation = reader[begin_translate+1:header_index]
	translation_dict = {}
	id_list = []
	for row in translation:
		trans = row[0].strip()
		tax_id = trans.split(' ')[0]
		tax_name = trans.split(' ')[1]
		id_list.append(tax_id)
		translation_dict.setdefault(tax_id,tax_name)
	
	return translation_dict,header_index,reader


def get_node_age_for_clade(clade,raw_tree,translation_dict):
	topology=[]
	for element in raw_tree:
		element = re.sub("^tree.+= ","",element)
		element = re.sub("\[.*?\]","",element)
		topology.append(element)
	topology = ','.join(topology)
	for element in translation_dict:
		# this regex says to replace any occurrence of the id in the tree, if preceded by any non-word character except for ':' and '-'. This is to ensure that not some numbers in the branchlength get replaced
		topology = re.sub(r'([^\w:-])%s(\b)' %element, r'\1%s\2' %translation_dict[element],topology, count=1)
	topology = re.sub(",\)",")",topology)
	topology = re.sub("\)\(","),(",topology)
	clean_topology = []
	tree_units = topology.split(':')
	for element in tree_units:
		element = re.sub("^[0-9-.E]+","",element,count=0)
		clean_topology.append(element)
	clean_topology = ''.join(clean_topology)

	final_match = ''
	#print clade
	#print clean_topology
	if clade in clean_topology:
		taxon_list = []
		for taxon in clade.split(','):
			taxon = re.sub('[()]','',taxon)
			taxon_list.append(taxon)
		taxa_w_regex = '.*'.join(taxon_list)
		search_string = '%s.*\)' %taxa_w_regex
		search_string = search_string+':[0-9-.E]+\)'       #the '\)' was added to the search pattern, it also worked without for the stacey results
		match = re.findall(search_string, topology)

		if len(match) == 1:
			match = match[0]
			#the number of closing parantheses for encompassing the clade is always 'the number of taxa'-1
			count = len(taxon_list)-1
			search_string2 = '.*?\)'*count
			final_match = re.findall(search_string2,match)
			final_match=final_match[0]
			#print final_match
		
		# retrieve all brlens information for all branches within the extracted clade
		branchlengths = re.findall(":[0-9.E-]+", final_match)
		
		# in order to calculate the node height of a clade, we need to weigh the internal brnaches by multiplying 
		#their values with a factor that corresponds to the number of terminals that are child to this branch. 
		#After weighing all brnaches in this matter, we can simply add all values together and divide it by the 
		#number of taxa, which will give us the brnachlength.

		#balance the parantheses in the match
		para_open = final_match.count('(')
		para_close = final_match.count(')')
		difference = para_close-para_open
		final_match = '('*difference + final_match

		#now go through the match starting from the back and find the enclosed clade
		final_brlengths = []
		#print branchlengths
		for brlens in branchlengths:
			brlen_match = re.findall('.*?%s'%brlens,topology)
			brlen_match = brlen_match[0]
			brlen_match = brlen_match.strip(brlens)
			reverse = brlen_match[::-1]
			counter = 0
			result = []
			for char in reverse:
				result.append(char)
				if char == ')':
					counter = counter+1
				if char == '(':
					counter = counter-1
				if counter == 0:
					#print 'reached'
					break
			factor = ''
			# if the result is only a single character it means that the brlens belonged to a tip in the tree and doesn't need to be weighed
			if len(result) == 1:
				factor = 1
			# in all other cases, determine how many taxa are enclosed by the clade that is defined by this branch
			else:
				final_result = ''.join(result)
				final_result = final_result[::-1]
				taxa_counter = 0
				for taxon in taxon_list:
					if taxon in final_result:
						taxa_counter += 1
				factor = taxa_counter
			final_brlengths.extend([brlens] * factor)

		converted_brlength = []
		for branchlenght in final_brlengths:
			converted_brlength.append(float(branchlenght.strip(':')))
		# now divide the sum of the final brlens dict by the number of taxa 
		node_height = sum(converted_brlength)/len(taxon_list)
		processed = 'true'
	else:
		node_height = 0
		processed = 'false'
	return node_height,processed



#######################################################################################################



translation_dict,header_index,reader = get_metadata(treefile)

body = reader[header_index+1:]
total_trees = len(body)-2
burnin_trees = int(total_trees*float(burnin))
print 'Skipping first',burnin_trees,'trees as burn-in'

if type(process_trees) == int:
	new_body = body[burnin_trees+1:burnin_trees+1+process_trees]
else:
	new_body = body[burnin_trees+1:]
print len(new_body), "trees are being analyzed"


nodeheight_dict = {}

for row in new_body:
	if len(row) > 0 and row[0] == 'End;':
		continue
	for clade in clade_list:
		node_height,processed = get_node_age_for_clade(clade,row,translation_dict)
		if processed == 'true':
			nodeheight_dict.setdefault(clade,[])
			nodeheight_dict[clade].append(node_height)

#print nodeheight_dict

for clade in nodeheight_dict:
	print 'The mean node heigt of clade %s is:' %clade,float(sum(nodeheight_dict[clade])/len(nodeheight_dict[clade])),'. Clade found in', len(nodeheight_dict[clade]), 'trees'

DE_output = open("./DE_node_depths.txt", "wb")
DE_output_log=csv.writer(DE_output, delimiter='\n')
DE_output_log.writerow(nodeheight_dict[clade_list[0]])

YZ_output = open("./YZ_node_depths.txt", "wb")
YZ_output_log=csv.writer(YZ_output, delimiter='\n')
YZ_output_log.writerow(nodeheight_dict[clade_list[1]])

XYZ_output = open("./XYZ_node_depths.txt", "wb")
XYZ_output_log=csv.writer(XYZ_output, delimiter='\n')
XYZ_output_log.writerow(nodeheight_dict[clade_list[2]])

DEXYZ_output = open("./DEXYZ_node_depths.txt", "wb")
DEXYZ_output_log=csv.writer(DEXYZ_output, delimiter='\n')
DEXYZ_output_log.writerow(nodeheight_dict[clade_list[3]])


































def get_clean_topology(raw_tree):
	topology = []
	for element in raw_tree:
		for id in id_list:
			element = re.sub("%s:[0-9.E-]+" %id,id,element)
		element = re.sub(":[0-9-.E]+","",element,count=0)			
		element = re.sub("^tree.+= ","",element)
		element = re.sub("\[.*?\]","",element)

		# this part is only relevant for more complex beast files with all sorts of elements in the tree file (such as substitution rates etc.)
		if '=' in element:
			if ')' in element:
				topology.append(')')
			if '(' in element:
				element = re.sub("^.+\(","(",element)
				topology.append(element)
		else:
			topology.append(element)

	topology = ','.join(topology)
	for element in translation_dict:
		# this regex says to replace any occurrence of the id in the tree, if preceded by any non-word character except for ':' and '-'. This is to ensure that not some numbers in the branchlength get replaced
		topology = re.sub(r'([^\w:-])%s(\b)' %element, r'\1%s\2' %translation_dict[element],topology, count=1)
	topology = re.sub(",\)",")",topology)
	topology = re.sub("\)\(","),(",topology)
	return(topology)


def get_branchlengths_topology(raw_tree):
	topology = []
	for element in row:
		element = re.sub("^tree.+= ","",element)
		element = re.sub("\[.*?\]","",element)

		# this part is only relevant for more complex beast files with all sorts of elements in the tree file (such as substitution rates etc.)
		if '=' in element:
			if ')' in element:
				topology.append(')')
			if '(' in element:
				element = re.sub("^.+\(","(",element)
				topology.append(element)
		else:
			topology.append(element)

	topology = ','.join(topology)
	for element in translation_dict:
		# this regex says to replace any occurrence of the id in the tree, if preceded by any non-word character except for ':' and '-'. This is to ensure that not some numbers in the branchlength get replaced
		topology = re.sub(r'([^\w:-])%s(\b)' %element, r'\1%s\2' %translation_dict[element],topology, count=1)
	topology = re.sub(",\)",")",topology)
	topology = re.sub("\)\(","),(",topology)
	return(topology)


## run through the tree file and find the most common topology (all topologies stored in dict)
#with open(treefile, 'r') as f:
#	reader = csv.reader(f, delimiter=',')
#	reader = list(reader)
#	header_index = ''
#	for row in reader:
#		if len(row) > 0 and row[0].startswith('tree STATE_0'):
#			header_index = int(reader.index(row)-1)
#			print 'skipping first', header_index, 'lines (header)'
#	
#	info = reader[0:header_index]
#	for row in info:
#		if len(row) > 0 and row[0].strip() == 'Translate':
#			begin_translate = info.index(row)
#	
#	translation = reader[begin_translate+1:header_index]
#	translation_dict = {}
#	id_list = []
#	for row in translation:
#		trans = row[0].strip()
#		tax_id = trans.split(' ')[0]
#		tax_name = trans.split(' ')[1]
#		id_list.append(tax_id)
#		translation_dict.setdefault(tax_id,tax_name)
#	#print translation_dict
#	
#	body = reader[header_index+1:]
#	total_trees = len(body)-2
#	burnin_trees = int(total_trees*float(burnin))
#	print 'Skipping first',burnin_trees,'trees as burn-in'
#
#	if process_trees != 'all':
#		new_body = body[burnin_trees+1:burnin_trees+1+process_trees]
#	elif process_trees == 'all':
#		new_body = body[burnin_trees+1:]
#	print len(new_body), "trees are being analyzed"
#	nodedepth_dict = {}
#
#
#	topology_dict = {}
#	for row in new_body:
#		if len(row) > 0 and row[0] == 'End;':
#			continue
#	
#		topology = get_clean_topology(row)
#
#		#make a dictionary of topologies to only choose the most common topology and skip all other lines
#		topology_dict.setdefault(topology,1)
#		if topology in topology_dict:
#			topology_dict[topology] = topology_dict[topology]+1
#
#sorted_topology_dict = sorted(topology_dict.items(), key=itemgetter(1))
##print sorted_topology_dict
#for topo in sorted_topology_dict:
#	if '((d1,d2),(e1,e2))' in topo[0]:
#		print topo
#most_common_topology = sorted_topology_dict[-1][0]
#
#
## now iterate through the file again, only selecting those trees representing the most common topology
#with open(treefile, 'r') as f:
#	reader = csv.reader(f, delimiter=',')
#	reader = list(reader)
#	header_index = ''
#	for row in reader:
#		if len(row) > 0 and row[0].startswith('tree STATE_0'):
#			header_index = int(reader.index(row)-1)
#	
#	info = reader[0:header_index]
#	for row in info:
#		if len(row) > 0 and row[0].strip() == 'Translate':
#			begin_translate = info.index(row)
#	
#	translation = reader[begin_translate+1:header_index]
#	translation_dict = {}
#	id_list = []
#	for row in translation:
#		trans = row[0].strip()
#		tax_id = trans.split(' ')[0]
#		tax_name = trans.split(' ')[1]
#		id_list.append(tax_id)
#		translation_dict.setdefault(tax_id,tax_name)
#	#print translation_dict
#	
#	body = reader[header_index+1:]
#	total_trees = len(body)-2
#	burnin_trees = int(total_trees*float(burnin))
#	#print 'Skipping first',burnin_trees,'trees as burn-in'
#	new_body = body[burnin_trees+1:]
#	#print len(new_body), "trees are being analyzed"
#	nodedepth_dict = {}
#
#
#	count = 0
#	for row in new_body:
#		if len(row) > 0 and row[0] == 'End;':
#			continue
#	
#		topology = get_clean_topology(row)
#		if topology == most_common_topology:
#
#			brlens_tree = get_branchlengths_topology(row)
#			branchlengths = re.findall(":[0-9.E-]+", brlens_tree)
#
#			if snapp != 1:
#				node_DE = float(branchlengths[0].strip(':'))+float(branchlengths[2].strip(':'))
#				node_YZ = float(branchlengths[10].strip(':'))
#				node_XYZ = float(branchlengths[12].strip(':'))+float(branchlengths[14].strip(':'))
#				root = float(branchlengths[12].strip(':'))+float(branchlengths[14].strip(':'))+float(branchlengths[15].strip(':'))
#			elif snapp == 1:
#				node_DE = float(branchlengths[9].strip(':'))+float(branchlengths[11].strip(':'))
#				node_YZ = float(branchlengths[3].strip(':'))
#				node_XYZ = float(branchlengths[5].strip(':'))+float(branchlengths[7].strip(':'))
#				root = float(branchlengths[5].strip(':'))+float(branchlengths[7].strip(':'))+float(branchlengths[8].strip(':'))			
#
#			nodedepth_dict.setdefault('DE',[])
#			nodedepth_dict['DE'].append(node_DE)
#			nodedepth_dict.setdefault('YZ',[])
#			nodedepth_dict['YZ'].append(node_YZ)
#			nodedepth_dict.setdefault('XYZ',[])
#			nodedepth_dict['XYZ'].append(node_XYZ)
#			nodedepth_dict.setdefault('DEXYZ',[])
#			nodedepth_dict['DEXYZ'].append(root)
#
#		else:
#			continue
#		#making sure that the correct number of trees is being processed
#		count += 1
#
#	print count, "trees with the most prominent topology are being summarized"
#
#	print 'The mean node heigt of clade DE is:',float(sum(nodedepth_dict['DE'])/len(nodedepth_dict['DE']))
#	print 'The mean node heigt of clade YZ is:',float(sum(nodedepth_dict['YZ'])/len(nodedepth_dict['YZ']))
#	print 'The mean node heigt of clade XYZ is:',float(sum(nodedepth_dict['XYZ'])/len(nodedepth_dict['XYZ']))
#	print 'The mean node heigt of clade DEXYZ is:',float(sum(nodedepth_dict['DEXYZ'])/len(nodedepth_dict['DEXYZ']))
#
#	DE_output = open("./DE_node_depths.txt", "wb")
#	DE_output_log=csv.writer(DE_output, delimiter='\n')
#	DE_output_log.writerow(nodedepth_dict['DE'])
#
#	YZ_output = open("./YZ_node_depths.txt", "wb")
#	YZ_output_log=csv.writer(YZ_output, delimiter='\n')
#	YZ_output_log.writerow(nodedepth_dict['YZ'])
#
#	XYZ_output = open("./XYZ_node_depths.txt", "wb")
#	XYZ_output_log=csv.writer(XYZ_output, delimiter='\n')
#	XYZ_output_log.writerow(nodedepth_dict['XYZ'])
#
#	DEXYZ_output = open("./DEXYZ_node_depths.txt", "wb")
#	DEXYZ_output_log=csv.writer(DEXYZ_output, delimiter='\n')
#	DEXYZ_output_log.writerow(nodedepth_dict['DEXYZ'])


		#nodedepth_dict.setdefault('DEXYZ',[])
		#nodedepth_dict['DEXYZ'].append(float(nodedepth[0].strip(':')))

	#print nodedepth_dict


#			if ']' in element or ':' in element:
#				#element = re.sub("\:.*[,)]",":",element)
#				print element
			
#			element = re.sub("[a-z]*[0-9]*\.rate=[0-9]*\.[0-9]*","",element)
#			element = re.sub(":\[&[a-z]*[0-9]*\.rate=.*","",element)

