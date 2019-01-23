import igraph
from igraph import *


#generate networks dictinary
def generate_networks():
    network_directory = "for_network_draw/cc_n100_65_14_d0.000885492883831.txt"
    network_files = open(network_directory,"r")

    networks={}
    lines = network_files.readlines()
    i=1
    for line in lines:
        nodes = {}
        network = "Network "+str(i)
        for gene in line.strip().split(" "):
            nodes[gene]=[]
        networks[network]=nodes
        i+=1
    return networks




#generate network edges
def generate_network_edges():
    edges={}
    edges_directory = "for_network_draw/hint_edge_file.txt"
    edges_file = open(edges_directory,"r")
    lines = edges_file.readlines()
    for line in lines:
        line =line.strip().split(" ")
        if line[0] in edges:
            edges[line[0]].append(line[1])
        else:
            edges[line[0]]=[]
            edges[line[0]].append(line[1])

        if line[1] in edges:
            edges[line[1]].append(line[0])
        else:
            edges[line[1]]=[]
            edges[line[1]].append(line[0])
    return edges



#generate gene_to_id, id_to_gene
def gen_2_id_and_id_2_gene():
    gene_name_to_id={}
    id_to_gene_name={}
    genes_directory = "for_network_draw/hint_index_file.txt"
    genes_file = open(genes_directory,"r")
    lines = genes_file.readlines()
    for line in lines:
        line =line.strip().split(" ")

        #print(line)
        gene_name_to_id[line[1].split("\t")[0]] = line[0]
        id_to_gene_name[line[0]] = line[1].split("\t")[0]

    return gene_name_to_id,id_to_gene_name

def generate_nodes_from_file(networks):
    net_nodes={}
    for network in networks:
        net_nodes[network] = []
        #print(networks[network])
        for gene in networks[network]:

            #print(net)
            gene_id = gene_equivqlence[gene]

            #print(gene_id)
            print(gene_id)
            gene_edges = edges[gene_id]
            print(gene_edges)
            #print("gene edges   ",gene_edges)
            #print(networks[network])

            for gene_ in gene_edges:

                print("iyaaa       ",gene_,"     gene_equivqlence_reversse[gene_]:    ",gene_equivqlence_reversse[gene_]   )

                #print(gene_)
                #print("here   ",gene_equivqlence_reversse[gene])
                #print(gene_equivqlence_reversse[gene] in net)
                #print (gene_equivqlence_reversse[gene_])
                #print(networks[network][gene])

                if (gene_equivqlence_reversse[gene_] in networks[network]):
                    if (str(gene_equivqlence_reversse[gene_]) not in networks[network][gene]):
                        networks[network][gene].append(gene_equivqlence_reversse[gene_])
                    if gene_equivqlence_reversse[gene_] not in net_nodes[network]:
                        net_nodes[network].append(gene_equivqlence_reversse[gene_])


            #print(networks[network])
        print("")
        print("")
        print("")
        print(networks[network])

    print(net_nodes)
    return net_nodes
#preparing data for igraph
def generate_igraph_tuples(networks,num_of_vertices):
    net_tuples_for_igraph={}
    last_vertex=num_of_vertices+1
    for net in networks:
        net_tuples_for_igraph[net]=[]
        for gene in networks[net]:
            for edge in networks[net][gene]:
                if net_nodes[net].index(gene) <= num_of_vertices and net_nodes[net].index(edge)<=num_of_vertices:
                    if ((net_nodes[net].index(edge),(net_nodes[net].index(gene))) not in net_tuples_for_igraph[net]):
                        net_tuples_for_igraph[net].append((net_nodes[net].index(gene),(net_nodes[net].index(edge))))


                elif net_nodes[net].index(gene) > num_of_vertices and net_nodes[net].index(edge)<=num_of_vertices:
                    if ((net_nodes[net].index(edge),last_vertex) not in net_tuples_for_igraph[net]) and (last_vertex,net_nodes[net].index(edge)) not in net_tuples_for_igraph[net]:
                        net_tuples_for_igraph[net].append((last_vertex,net_nodes[net].index(edge)))

                elif net_nodes[net].index(gene) <= num_of_vertices and net_nodes[net].index(edge)> num_of_vertices:
                    if ((net_nodes[net].index(gene),last_vertex) not in net_tuples_for_igraph[net]) and (last_vertex,net_nodes[net].index(gene)) not in net_tuples_for_igraph[net]:
                        net_tuples_for_igraph[net].append((last_vertex,net_nodes[net].index(gene)))

    return net_tuples_for_igraph


num_of_vertices=9
network_num="Network 12"
networks= generate_networks()
edges = generate_network_edges()
gene_equivqlence,gene_equivqlence_reversse = gen_2_id_and_id_2_gene()
net_nodes = generate_nodes_from_file(networks)
net_tuples_for_igraph = generate_igraph_tuples(networks,num_of_vertices)




g = Graph(net_tuples_for_igraph[network_num], directed=None,vertex_attrs={'y': net_nodes[network_num]})
net_nodes[network_num][num_of_vertices+1] = "Others ["+ str(len(net_nodes[network_num])-num_of_vertices) +"]"
print(net_nodes[network_num])
g.vs["label"] = net_nodes[network_num]
visual_style = {}
visual_style["vertex_size"] = 50
visual_style["vertex_color"] = "green"
visual_style["layout"] = "circle"
visual_style["vertex_shape"] = "rectangle"
visual_style["bbox"] = (800, 800)
visual_style["margin"] = 50
plot(g, **visual_style)
#layout = g.layout("kk")
#plot(g, vertex_shape="rectangle", vertex_size=50,vertex_color="green", layout = "circle")
