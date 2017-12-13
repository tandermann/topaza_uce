import os
import io
import re
import csv
from Bio import Phylo


#treefile = raw_input("Give path to tree file: ")
treefile = "/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/simulated/allele_alignments/rescaled_species.trees"
subfolder = treefile.split('/')[-3]+'_'+treefile.split('/')[-2]
out_dir_raw = "/Users/tobias/Desktop/node_depth_distribution/allele_alignments"
out_dir = os.path.join(out_dir_raw,subfolder)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
burnin = 0.1
process_trees = 'all'
clade_list = ['d,e','y,z','x,y,z','e,d,x,y,z']
#######################################################################################################


def get_metadata(treefile):
    reader = csv.reader(open(treefile, 'r'), delimiter=',')
    reader = list(reader)
    header_index = ''
    for row in reader:
        if len(row) > 0 and row[0].startswith('tree STATE_0'):
            header_index = int(reader.index(row)-1)
            print('Skipping first', header_index, 'lines (header)')
        
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


def replace_numbers_with_taxon_names(raw_tree,translation_dict):
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
    return topology


def get_tip_dict(tree_obj,clades):
    brlens_dict = tree_obj.depths(0)
    # Get all terminal clades (=tips) belonging to each species/taxon
    clade_dict = {}
    for clade in brlens_dict.keys():
        reduced_name = None
        if clade.name:
            # Here the taxon name is simply the first letter of every tip name, so reduce the tip name to it's first letter
            reduced_name = re.sub('^(.).*','\\1',clade.name)
        if reduced_name:
            match = None
            for target_clade in clades:
                if reduced_name in target_clade.split(','):
                    match = target_clade
                if match:
                    clade_dict.setdefault(match,[])
                    if not clade in clade_dict[match]:
                        clade_dict[match].append(clade)
    return clade_dict


def mrca_node_depth_for_clades(tree_object,clade_dict):
    # Get the common ancestor for the tips belonging to each species/taxon
    #root_depth_dict = {}
    taxon_crown_height_dict = {}
    for taxon in clade_dict.keys():
        #root_depth_dict.setdefault(taxon,[])
        taxon_crown_height_dict.setdefault(taxon,[])
        common_ancestor = tree_object.common_ancestor(clade_dict[taxon])
        # This function retrieves all clades that are between the common ancesotr root and the tips. We want the brlens of all of these branches
        #involved_clades = tree_object.trace(clade_dict[taxon][0], common_ancestor)
        #root_length = common_ancestor.branch_length
        for tip in clade_dict[taxon]:
            # for every tip we want to store the length of the complete path, starting at the root
            #tip_total_length = []
            taxon_crown_age = []
            #tip_total_length.append(root_length)
            #tip_name = tip.name
            parent_clades = common_ancestor.get_path(target=tip)
            for clade in parent_clades:
                #tip_total_length.append(clade.branch_length)
                taxon_crown_age.append(clade.branch_length)
            #tip_total_length = sum(tip_total_length)
            #root_depth_dict[taxon].append(tip_total_length)
            taxon_crown_age = sum(taxon_crown_age)
            taxon_crown_height_dict[taxon].append(taxon_crown_age)
    return taxon_crown_height_dict




#######################################################################################################

# removing header and burnin form trees file:
translation_dict,header_index,reader = get_metadata(treefile)
body = reader[header_index+1:]
total_trees = len(body)-2
burnin_trees = int(total_trees*float(burnin))
print('Skipping first',burnin_trees,'trees as burn-in')
if type(process_trees) == int:
	new_body = body[burnin_trees+1:burnin_trees+1+process_trees]
else:
	new_body = body[burnin_trees+1:]
print(len(new_body), "trees of in total", len(body),"are being analyzed")


nodeheight_dict = {}
# iterate through trees
for row in new_body:
    topology = replace_numbers_with_taxon_names(row,translation_dict)
    tree_object = Phylo.read(io.StringIO(topology), "newick")
    #Phylo.draw(tree_object)
    #print(tree_object)
    tip_dict = get_tip_dict(tree_object,clade_list)
    node_depth_dict = mrca_node_depth_for_clades(tree_object,tip_dict)    
    for clade in node_depth_dict.keys():
        # get one of the values (all identical) for the node height for each clade 
        node_height = node_depth_dict[clade][0]
        nodeheight_dict.setdefault(clade,[])
        nodeheight_dict[clade].append(node_height)


for clade in nodeheight_dict:
    print('The mean node heigt of clade %s is:' %clade,float(sum(nodeheight_dict[clade])/len(nodeheight_dict[clade])))

DE_output = open(os.path.join(out_dir,"./DE_node_depths.txt"), "w")
DE_output_log=csv.writer(DE_output, delimiter='\n')
DE_output_log.writerow(nodeheight_dict[clade_list[0]])

YZ_output = open(os.path.join(out_dir,"./YZ_node_depths.txt"), "w")
YZ_output_log=csv.writer(YZ_output, delimiter='\n')
YZ_output_log.writerow(nodeheight_dict[clade_list[1]])

XYZ_output = open(os.path.join(out_dir,"./XYZ_node_depths.txt"), "w")
XYZ_output_log=csv.writer(XYZ_output, delimiter='\n')
XYZ_output_log.writerow(nodeheight_dict[clade_list[2]])

DEXYZ_output = open(os.path.join(out_dir,"./DEXYZ_node_depths.txt"), "w")
DEXYZ_output_log=csv.writer(DEXYZ_output, delimiter='\n')
DEXYZ_output_log.writerow(nodeheight_dict[clade_list[3]])


