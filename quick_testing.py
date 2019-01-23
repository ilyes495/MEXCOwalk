with open('../hint/out/connected_components_isolarge/memcover_v3/cc_n1100_memcover_hint_0.03_k3_genes_countrealN1100.txt', 'r') as f:
    lines = f.readlines()
    genes =[]
    for l in lines:
        genes.extend(l.rstrip().split('\t'))
    print(genes[:20])
    print(len(genes), len(set(genes)))
