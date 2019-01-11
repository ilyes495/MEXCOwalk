import glob
import sys
import math


godepthfile = open("../analysis_files/go.obo.dag.txt.level5.txt")
goassocfile = open("../analysis_files/hs_goassociations17")


# two input files
hotnet_file = sys.argv[1]
our_file = sys.argv[2]
fhout = open(sys.argv[3],'a')
files = [hotnet_file, our_file ]
hotnet_file_core = hotnet_file[hotnet_file.rindex('/'):]
num_genes = hotnet_file_core.split('_')[1][1:] #******I automatically find k from hotnet2 file, this might need to be modified accordingly.****
print files


godepthdic = {}
for line in godepthfile:
	cols = line.split()
	godepthdic[cols[0]] = set()
	for ind in range(1, len(cols)):
		if cols[ind][0] == 'G':
			godepthdic[cols[0]].add(cols[ind])


togo = {}
for line in goassocfile:
	cols = line.split()
	curgoset = set()
	nameset = set()
	for name in cols:
		if name[0] == 'G' and name in godepthdic:
			curgoset = curgoset.union(godepthdic[name])
		else:
			if name[0] != 'G' or name[1] != 'O':
				nameset.add(name)

	if len(curgoset) > 0:
		for name in nameset:
			if name not in togo:
				togo[name] = set()
			togo[name] = togo[name].union(curgoset)


for i in range(len(files)):
  filename = files[i]
  file_sum = 0
  fhinput = open(filename)
  module_count = 0
  module_scores = []
  for line in fhinput:
      module_count += 1
      genes = line.split()

      print genes
      all_scores = []

      count_genes_with_go = 0
      union_go = set()
      dict_go_gene = {}
      for gene in genes:
          if gene in togo:
              count_genes_with_go += 1
              union_go = union_go.union(togo[gene])
              for go_term in togo[gene]:
                   if go_term not in dict_go_gene:
                      dict_go_gene[go_term] = []
                   dict_go_gene[go_term].append(gene)

      num_go_terms = len(union_go)           # d in the formula
      if num_go_terms > 1: # if it's 1 the norm_term is a problem because of division by 0
          norm_term = - 1.0 / math.log(num_go_terms)
          sum_ne = 0
          for go_term in union_go:
             #print ' go term ', len(dict_go_gene[go_term])
             frac_go_term = len(dict_go_gene[go_term]) / float(count_genes_with_go)
             if frac_go_term == 0:
                sum_ne += 0
             else:
                sum_ne += frac_go_term * math.log(frac_go_term)

          module_score = sum_ne * norm_term
          module_scores.append(module_score)
      #print 'len genes ', len(genes), ' genes with go ', count_genes_with_go, ' num go terms ', num_go_terms, ' sum ne ', sum_ne, ' norm term ', norm_term, ' module score ', module_score

  final_score = sum(module_scores) / len(module_scores)
  print >>fhout, num_genes
  if i == 0:
       print >>fhout, round(final_score, 4), "our"
  else:
       print >>fhout, round(final_score, 4), "hotnet"
  print >>fhout, '****'


fhout.close()
