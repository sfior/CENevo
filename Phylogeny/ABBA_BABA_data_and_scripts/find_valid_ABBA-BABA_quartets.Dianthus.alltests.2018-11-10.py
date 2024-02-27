from ete2 import Tree
import ete2
import random


def make_nonredundant_pairs (samples):
	
	iterlist1 = samples
	iterlist2 = iterlist1
	
	combinations = []
	for a in iterlist1:
		for b in iterlist2:
			if a != b:
				comb = sorted([a, b])
				combinations.append(comb)
	
	unique_combs = [list(x) for x in set(tuple(x) for x in combinations)]

	return unique_combs

def get_tip_names ( t_node ):
	
	leaves = t_node.get_leaves(is_leaf_fn=None)
	tips = []
	for l in leaves:
		tips.append( l.name)
	
	return tips

def cull_down_quartets ( good_quartets, max_n_contrasts ):
	
	# this function respects the directionality of the comparison:
	# outgr, A, B, C | outgr, B, A, D may both be valid quartets on the overall topology,
	# and both test migration between A and B,
	# but differ in their direction of geneflow (if any)
	
	# this works only if all quartets have the same outgroup
	outgr = good_quartets[0][0]
	
	contrasts_dict = {}
	for q in good_quartets:
		leaf3, leaf2, leaf1 = q[1], q[2], q[3]
		try:
			contrasts_dict[ "$".join(sorted( [leaf3, leaf2] )) ].append( q )
		except KeyError:
			contrasts_dict[ "$".join(sorted( [leaf3, leaf2] )) ] = [ q ]
	
	# the downsampling:
	culled_quartets = []
	for pair in contrasts_dict.keys():
		try:
			accepted_contrasts = random.sample(contrasts_dict[pair], max_n_contrasts)
		except ValueError:
			accepted_contrasts = contrasts_dict[pair]
		for a in accepted_contrasts:
			culled_quartets.append( a )
	
	print """after downsampling to max. {0} contrasts per species pair,
retained {1} ABBA-BABA quartets ( {2} % )""".format( max_n_contrasts , len( culled_quartets ) , float(len(culled_quartets))/float(len( good_quartets ))*100 )
		
	return culled_quartets

def make_quartets_treewalking (t, outgr):
	"""	
		# idea: starting from outgr, get next node, get all 
		for internal_node in allnodes:
			if len(leaves_from_node) >= 3:
				left_leaves =
				right_leaves =
				all_left_pairs = nonredundant_pairs
				all_right_pairs = nonredundant_pairs
				left_triplets = all_left_pairs * right_leaves
				right_triplets = all_right_pairs * left_leaves
			
				quartets.append( [  ] ) 
	"""
	
	all_valid_quartets = []
	for node in t.get_descendants(): # iterate over all internal (and terminal = leaf) nodes
		if len( node.get_descendants() ) >= 3: # reject leaves and nodes that have too few leaves on them
			# distinguish the two descendant subtrees:
			left_subtree = node.children[0]
 			right_subtree = node.children[1]
			
			left_tips = get_tip_names( left_subtree )
			right_tips = get_tip_names( right_subtree )
			
			# get all possible leaf2, leaf1 pairings:
			all_left_pairs = make_nonredundant_pairs ( left_tips )
			all_right_pairs = make_nonredundant_pairs ( right_tips )
			
			# get the triplets:
			left_triplets = []
			for l in all_left_pairs:
				for r in right_tips:
					all_valid_quartets.append( [outgr, r, l[0], l[1] ] )
			
			right_triplets = []
			for r in all_right_pairs:
				for l in left_tips:
					all_valid_quartets.append( [outgr, l, r[0], r[1] ] )
			
# 			print left_subtree
# 			print left_triplets
# 			print right_subtree
# 			print right_triplets
# 			print "-----------------------------"
	
	bad_samples= ["sample_29_D_deltoides","sample_50_D_longicaulis","sample_4_D_arenarius","sample_44_D_glacialis","sample_22_D_carthusianorum",
"sample_25_D_carthusianorum",
"sample_27_D_carthusianorum",
"sample_26_D_carthusianorum",
"sample_43_D_glacialis",
"sample_86_D_viscidus",
"sample_76_D_ssp.",
"sample_13_D_broteri",
"sample_16_D_carthusianorum",
"sample_12_D_brachycalyx",
"sample_52_D_longicaulis",
"sample_17_D_carthusianorum",
"sample_11_D_brachycalyx",
"sample_Ds_POP3",
"sample_48_D_hungaricus",
"sample_78_D_strictus",
"sample_49_D_juniperinus_subsp_bauhinorum",
"sample_Dc_POP6",
"sample_Dc_POP4",
"sample_Dc_POP3",
"sample_Dc_POP5",
"sample_Dc_POP2",
"sample_Dc_POP1"]
	
	selected_quartets = []
	for q in all_valid_quartets:
		good = True
		armerium_count = 0
		for sp in q:
			if "armeria" in sp or "deltoides" in sp or "viscidus" in sp or "D_ssp" in sp:
				armerium_count +=1
		if armerium_count > 1:
			good = False 	
		for sp in q:
			for bad_sample in bad_samples:
				if bad_sample in q:
					good = False				
		if good:
			selected_quartets.append(q)
	
	
	print "found {0} valid ABBA-BABA quartets".format( len( selected_quartets ) )
	
	return selected_quartets

def sort_and_output ( quartets ):
	
	## sort the outcome and write to file::
	quart_dict = {}
	for i in culled_quartets:
		quart_dict["".join(i)] = i
	
#	outlines = ["\t".join(["outgr", "P3", "P2", "P1"]) ]
	outlines = []
	for i in sorted(quart_dict.keys()):
		quartet = quart_dict[ i ]
		quartet.reverse()
		outlines.append( "\t".join( quartet ))
	
	with open("ABBA-BABA_quartets.Dianthus.alltests.2018-11-10.txt", "a") as OUTFILE:
		OUTFILE.write("\n".join(outlines) + "\n")

			
############# MAIN

treefile = "Dianthus_tree.descr_names.tre"
#outgr = "sample_9_D_armeria" # higher coverage
#outgr = "sample_34_D_deltoides" # highest coverage
outgr = "sample_87_D_viscidus" # highest coverage

max_n_contrasts = 10



t = Tree( treefile )
t.set_outgroup( outgr )

all_valid_quartets = make_quartets_treewalking (t, outgr)

culled_quartets = cull_down_quartets ( all_valid_quartets, max_n_contrasts )

sort_and_output ( culled_quartets )

print "Done!"
