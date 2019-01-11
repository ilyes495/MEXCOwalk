import numpy as np

GENE_ID_OFFSET = 1  # the genes id starts with 1

network_name = "hint"     # ["hint", "irefindex", "multinet"]
network_beta = 0.4        # [0.4, 0.45, 0.5]

def load_gene_vs_patient_data():
    data = {}
    with open("../data/genes_vs_patients_indices_gen_paran.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1:]
    return data
    
def load_gene_vs_patient_data_ready():
    data = {}
    with open("../data/pan12gene2freq.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = float(line[1])
    return data

def load_patients_to_indices():
    data = {}
    with open("../data/patients_to_indices_gen.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1]
    return data

def load_src_cnas_data():
    data = {}
    with open("../data/src_cnas.tsv") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1:]
    return data

def load_src_snvs_data():
    data = {}
    with open("../data/src_snvs.tsv") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1:]
    return data

def load_id_to_gene():
    filename = ""
    if network_name == "hint":
        filename = "../hint/data/hint_index_file.txt"

    id_to_gene = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[int(line[0]) - GENE_ID_OFFSET] = line[1]
    return id_to_gene

def load_gene_list():
    filename = ""
    if network_name == "hint":
        filename = "../hint/data/hint_index_file.txt"

    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[1])
    return genes

def load_gene_to_id():
    filename = ""
    if network_name == "hint":
        filename = "../hint/data/hint_index_file.txt"

    id_to_gene = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[line[1]] = int(line[0]) - GENE_ID_OFFSET
    return id_to_gene

def load_unique_genes():
    filename = ""
    if network_name == "hint":
        filename = "../hint/data/hint_index_file.txt"

    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[1])
    return genes

def load_edge_list():
    filename = ""
    if network_name == "hint":
        filename = "../hint/data/hint_edge_file.txt"

    edge_list = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            edge_list.append([int(line[0]) - GENE_ID_OFFSET, int(line[1]) - GENE_ID_OFFSET])
    return edge_list

def write_matrix_file(m, filename):
    with open(filename, "w+") as f:
        for i in range(len(m)):
            for j in range(len(m[i])):
                f.write(str(m[i][j])+" ")
            f.write("\n")
    return

def load_matrix_file(filename):
    n = len(load_gene_list()) # number of genes
    m = np.zeros((n, n))
    with open(filename) as f:
        lines = f.readlines()
        i = 0
        for line in lines:
            line = line.split()
            for j in range(n):
                m[i][j] = line[j]
            i = i + 1
    return m

hotnet_filenames = [
"cc_n100_k3_d0.000395306627809.txt",
"cc_n200_k3_d0.000317462997192.txt",
"cc_n300_k3_d0.000265161413542.txt",
"cc_n400_k3_d0.000237365057565.txt",
"cc_n500_k3_d0.000218844122647.txt",
"cc_n600_k3_d0.000197718399657.txt",
"cc_n700_k3_d0.000179260103238.txt",
"cc_n800_k3_d0.000168712691679.txt",
"cc_n900_k3_d0.00015689340613.txt",
"cc_n1000_k3_d0.000145393673689.txt",
"cc_n1100_k3_d0.000136391282202.txt",
"cc_n1200_k3_d0.000127520440545.txt",
"cc_n1300_k3_d0.000120673154974.txt",
"cc_n1400_k3_d0.000113156917812.txt",
"cc_n1500_k3_d0.000107632074618.txt",
"cc_n1600_k3_d0.000102426336924.txt",
"cc_n1700_k3_d9.86836417676e-05.txt",
"cc_n1800_k3_d9.50965371765e-05.txt",
"cc_n1900_k3_d9.13171011891e-05.txt",
"cc_n2000_k3_d8.7914150296e-05.txt",
"cc_n2100_k3_d8.50186314199e-05.txt",
"cc_n2200_k3_d8.1555942637e-05.txt",
"cc_n2300_k3_d7.75666930716e-05.txt",
"cc_n2400_k3_d7.35454213437e-05.txt",
"cc_n2500_k3_d6.95067334321e-05.txt",
"cc_n2600_k3_d6.62011987272e-05.txt",
"cc_n2700_k3_d6.30842875967e-05.txt",
"cc_n2800_k3_d6.06669640213e-05.txt",
"cc_n2900_k3_d5.84775749375e-05.txt",
"cc_n3000_k3_d5.56411985723e-05.txt",
"cc_n3100_k3_d5.34954488723e-05.txt",
"cc_n3200_k3_d5.14988020615e-05.txt",
"cc_n3300_k3_d4.99018354394e-05.txt",
"cc_n3400_k3_d4.84344260291e-05.txt",
"cc_n3500_k3_d4.58564718931e-05.txt",
"cc_n3600_k3_d4.37636627794e-05.txt",
"cc_n3700_k3_d4.2205206617e-05.txt",
"cc_n3800_k3_d4.05526780107e-05.txt",
"cc_n3900_k3_d3.88056787704e-05.txt",
"cc_n4000_k3_d3.67357705895e-05.txt",
"cc_n4100_k3_d3.52136077485e-05.txt",
"cc_n4200_k3_d3.37788734755e-05.txt",
"cc_n4300_k3_d3.20953438704e-05.txt",
"cc_n4400_k3_d3.06688026148e-05.txt",
"cc_n4500_k3_d2.9302895779e-05.txt",
"cc_n4600_k3_d2.78041001469e-05.txt",
"cc_n4700_k3_d2.63205426338e-05.txt",
"cc_n4800_k3_d2.50703923968e-05.txt",
"cc_n4900_k3_d2.40315156784e-05.txt",
"cc_n5000_k3_d2.28049555608e-05.txt",
"cc_n5100_k3_d2.1184001491e-05.txt",
"cc_n5200_k3_d1.99402094056e-05.txt",
"cc_n5300_k3_d1.89164060983e-05.txt"
]


our_filenames_a = [
"cc_n100_k3_d0.000397873369233.txt",
"cc_n200_k3_d0.000320233046807.txt",
"cc_n300_k3_d0.000269269823807.txt",
"cc_n400_k3_d0.000238263241126.txt",
"cc_n500_k3_d0.00022184060303.txt",
"cc_n600_k3_d0.000200119685724.txt",
"cc_n700_k3_d0.000181296272946.txt",
"cc_n800_k3_d0.000169816267025.txt",
"cc_n900_k3_d0.000158193643053.txt",
"cc_n1000_k3_d0.000146051721339.txt",
"cc_n1100_k3_d0.000136919060195.txt",
"cc_n1200_k3_d0.000128190863711.txt",
"cc_n1300_k3_d0.000121478002225.txt",
"cc_n1400_k3_d0.00011392343965.txt",
"cc_n1500_k3_d0.000109575413618.txt",
"cc_n1600_k3_d0.000103615790537.txt",
"cc_n1700_k3_d9.92244107456e-05.txt",
"cc_n1800_k3_d9.59137199981e-05.txt",
"cc_n1900_k3_d9.26522015213e-05.txt",
"cc_n2000_k3_d8.8905128953e-05.txt",
"cc_n2100_k3_d8.56298975074e-05.txt",
"cc_n2200_k3_d8.25458688383e-05.txt",
"cc_n2300_k3_d7.90911491358e-05.txt",
"cc_n2400_k3_d7.43042078918e-05.txt",
"cc_n2500_k3_d7.04736227734e-05.txt",
"cc_n2600_k3_d6.74217249923e-05.txt",
"cc_n2700_k3_d6.40524137086e-05.txt",
"cc_n2800_k3_d6.10335780964e-05.txt",
"cc_n2900_k3_d5.90547957546e-05.txt",
"cc_n3000_k3_d5.59559270818e-05.txt",
"cc_n3100_k3_d5.39183060972e-05.txt",
"cc_n3200_k3_d5.20231152826e-05.txt",
"cc_n3300_k3_d5.01665700038e-05.txt",
"cc_n3400_k3_d4.87182055096e-05.txt",
"cc_n3500_k3_d4.6534323777e-05.txt",
"cc_n3600_k3_d4.43472672821e-05.txt",
"cc_n3700_k3_d4.24898523009e-05.txt",
"cc_n3800_k3_d4.10325477433e-05.txt",
"cc_n3900_k3_d3.90709635207e-05.txt",
"cc_n4000_k3_d3.71449111584e-05.txt",
"cc_n4100_k3_d3.5543112911e-05.txt",
"cc_n4200_k3_d3.40779869591e-05.txt",
"cc_n4300_k3_d3.2293869704e-05.txt",
"cc_n4400_k3_d3.09145554733e-05.txt",
"cc_n4500_k3_d2.94280795388e-05.txt",
"cc_n4600_k3_d2.79280581362e-05.txt",
"cc_n4700_k3_d2.66061884309e-05.txt",
"cc_n4800_k3_d2.52816783187e-05.txt",
"cc_n4900_k3_d2.42749330633e-05.txt",
"cc_n5000_k3_d2.30710951099e-05.txt",
"cc_n5100_k3_d2.15250018476e-05.txt",
"cc_n5200_k3_d2.00081095676e-05.txt",
"cc_n5300_k3_d1.89581559746e-05.txt"
]

our_filenames_b = [
"cc_n100_k3_d0.000872929945757.txt",
"cc_n200_k3_d0.000596135951343.txt",
"cc_n300_k3_d0.000456207544277.txt",
"cc_n400_k3_d0.000384739789007.txt",
"cc_n500_k3_d0.000347869381092.txt",
"cc_n600_k3_d0.000308863525976.txt",
"cc_n700_k3_d0.00027624069105.txt",
"cc_n800_k3_d0.000245605987923.txt",
"cc_n900_k3_d0.000221556571159.txt",
"cc_n1000_k3_d0.000199347580934.txt",
"cc_n1100_k3_d0.000184281059152.txt",
"cc_n1200_k3_d0.000171576532535.txt",
"cc_n1300_k3_d0.000161368064118.txt",
"cc_n1400_k3_d0.000151950121201.txt",
"cc_n1500_k3_d0.000143150397691.txt",
"cc_n1600_k3_d0.000133356057692.txt",
"cc_n1700_k3_d0.000125255959129.txt",
"cc_n1800_k3_d0.000117869276268.txt",
"cc_n1900_k3_d0.000111411250882.txt",
"cc_n2000_k3_d0.000104956830649.txt",
"cc_n2100_k3_d9.9904653813e-05.txt",
"cc_n2200_k3_d9.44318916845e-05.txt",
"cc_n2300_k3_d8.98098900903e-05.txt",
"cc_n2400_k3_d8.45346374884e-05.txt",
"cc_n2500_k3_d8.0956055881e-05.txt",
"cc_n2600_k3_d7.73243348375e-05.txt",
"cc_n2700_k3_d7.37680581262e-05.txt",
"cc_n2800_k3_d7.04636118929e-05.txt",
"cc_n2900_k3_d6.61881430027e-05.txt",
"cc_n3000_k3_d6.34393294183e-05.txt",
"cc_n3100_k3_d6.05406338993e-05.txt",
"cc_n3200_k3_d5.76797001726e-05.txt",
"cc_n3300_k3_d5.42757762101e-05.txt",
"cc_n3400_k3_d5.19597742955e-05.txt"
]

our_filenames_c = [
"cc_n100_k3_d0.000478470396555.txt",
"cc_n200_k3_d0.000355251685457.txt",
"cc_n300_k3_d0.000287349093411.txt",
"cc_n400_k3_d0.000253797793454.txt",
"cc_n500_k3_d0.000233936934971.txt",
"cc_n600_k3_d0.000216279537899.txt",
"cc_n700_k3_d0.000194135569264.txt",
"cc_n800_k3_d0.000180957167701.txt",
"cc_n900_k3_d0.000169561961341.txt",
"cc_n1000_k3_d0.00015815144782.txt",
"cc_n1100_k3_d0.000147882326297.txt",
"cc_n1200_k3_d0.000137524359339.txt",
"cc_n1300_k3_d0.000129769116715.txt",
"cc_n1400_k3_d0.000124324231041.txt",
"cc_n1500_k3_d0.000117213808998.txt",
"cc_n1600_k3_d0.000111283457306.txt",
"cc_n1700_k3_d0.000105128146463.txt",
"cc_n1800_k3_d0.000100658674394.txt",
"cc_n1900_k3_d9.67577057271e-05.txt",
"cc_n2000_k3_d9.32051444471e-05.txt",
"cc_n2100_k3_d8.94539372944e-05.txt",
"cc_n2200_k3_d8.6541138098e-05.txt",
"cc_n2300_k3_d8.30167124818e-05.txt",
"cc_n2400_k3_d7.90655189548e-05.txt",
"cc_n2500_k3_d7.48228667755e-05.txt",
"cc_n2600_k3_d7.10560587115e-05.txt",
"cc_n2700_k3_d6.81715506164e-05.txt",
"cc_n2800_k3_d6.46993432822e-05.txt",
"cc_n2900_k3_d6.18684207371e-05.txt",
"cc_n3000_k3_d5.91981334182e-05.txt",
"cc_n3100_k3_d5.66975363628e-05.txt",
"cc_n3200_k3_d5.4061653936e-05.txt",
"cc_n3300_k3_d5.21685037107e-05.txt",
"cc_n3400_k3_d5.01354197824e-05.txt",
"cc_n3500_k3_d4.86687413078e-05.txt",
"cc_n3600_k3_d4.69308999751e-05.txt",
"cc_n3700_k3_d4.4508224912e-05.txt",
"cc_n3800_k3_d4.27246207788e-05.txt"
]

our_filenames_d = [
"cc_n100_k3_d0.00039407790271.txt",
"cc_n200_k3_d0.000317269860144.txt",
"cc_n300_k3_d0.00026455527812.txt",
"cc_n400_k3_d0.000235484645259.txt",
"cc_n500_k3_d0.00021799220487.txt",
"cc_n600_k3_d0.000198360527235.txt",
"cc_n700_k3_d0.000178738914549.txt",
"cc_n800_k3_d0.000168692754942.txt",
"cc_n900_k3_d0.000157011357198.txt",
"cc_n1000_k3_d0.000144921344854.txt",
"cc_n1100_k3_d0.000136355265014.txt",
"cc_n1200_k3_d0.000127654261361.txt",
"cc_n1300_k3_d0.00012075550688.txt",
"cc_n1400_k3_d0.000113184095001.txt",
"cc_n1500_k3_d0.000108143102835.txt",
"cc_n1600_k3_d0.000102595107524.txt",
"cc_n1700_k3_d9.87307146474e-05.txt",
"cc_n1800_k3_d9.50609859507e-05.txt",
"cc_n1900_k3_d9.17600014931e-05.txt",
"cc_n2000_k3_d8.78964186094e-05.txt",
"cc_n2100_k3_d8.48927843022e-05.txt",
"cc_n2200_k3_d8.15786828389e-05.txt",
"cc_n2300_k3_d7.78681700997e-05.txt",
"cc_n2400_k3_d7.35744499437e-05.txt",
"cc_n2500_k3_d6.95478231493e-05.txt",
"cc_n2600_k3_d6.63465486504e-05.txt",
"cc_n2700_k3_d6.31457252892e-05.txt",
"cc_n2800_k3_d6.06668418675e-05.txt",
"cc_n2900_k3_d5.86419995902e-05.txt",
"cc_n3000_k3_d5.57101606122e-05.txt",
"cc_n3100_k3_d5.35422961291e-05.txt",
"cc_n3200_k3_d5.16662504114e-05.txt",
"cc_n3300_k3_d4.98598757954e-05.txt",
"cc_n3400_k3_d4.83396390161e-05.txt",
"cc_n3500_k3_d4.58427530439e-05.txt",
"cc_n3600_k3_d4.37358298443e-05.txt",
"cc_n3700_k3_d4.21724121196e-05.txt",
"cc_n3800_k3_d4.05860828619e-05.txt"
]

our_filenames_e = [
"cc_n100_k3_d0.000666293563885.txt",
"cc_n200_k3_d0.000545533062112.txt",
"cc_n300_k3_d0.000426934503735.txt",
"cc_n400_k3_d0.000356339627396.txt",
"cc_n500_k3_d0.000311477937833.txt",
"cc_n600_k3_d0.0002826015996.txt",
"cc_n700_k3_d0.000260436271105.txt",
"cc_n800_k3_d0.000248852617401.txt",
"cc_n900_k3_d0.00022862830474.txt",
"cc_n1000_k3_d0.000210688502484.txt",
"cc_n1100_k3_d0.000196114183143.txt",
"cc_n1200_k3_d0.000184223022399.txt",
"cc_n1300_k3_d0.000175454907439.txt",
"cc_n1400_k3_d0.000165008629135.txt",
"cc_n1500_k3_d0.000151529080254.txt",
"cc_n1600_k3_d0.000144167281692.txt",
"cc_n1700_k3_d0.000132315670678.txt",
"cc_n1800_k3_d0.000125593158507.txt",
"cc_n1900_k3_d0.000118770986566.txt",
"cc_n2000_k3_d0.000113478978483.txt",
"cc_n2100_k3_d0.00010910891001.txt",
"cc_n2200_k3_d0.000103640200528.txt",
"cc_n2300_k3_d9.95078130514e-05.txt",
"cc_n2400_k3_d9.58795780387e-05.txt",
"cc_n2500_k3_d9.22286960234e-05.txt",
"cc_n2600_k3_d8.91698857336e-05.txt",
"cc_n2700_k3_d8.57990401145e-05.txt",
"cc_n2800_k3_d8.21097651247e-05.txt",
"cc_n2900_k3_d7.80619915824e-05.txt",
"cc_n3000_k3_d7.48439772688e-05.txt",
"cc_n3100_k3_d7.17127898109e-05.txt",
"cc_n3200_k3_d6.78859655367e-05.txt",
"cc_n3300_k3_d6.45618790627e-05.txt",
"cc_n3400_k3_d6.12098835134e-05.txt",
"cc_n3500_k3_d5.84154765598e-05.txt",
"cc_n3600_k3_d5.54399322697e-05.txt",
"cc_n3700_k3_d5.31576300012e-05.txt",
"cc_n3800_k3_d5.12120577246e-05.txt"
]

our_filenames_f = [
"cc_n100_k3_d0.000252640259556.txt",
"cc_n200_k3_d0.000194285590017.txt",
"cc_n300_k3_d0.000171559763584.txt",
"cc_n400_k3_d0.000145618799612.txt",
"cc_n500_k3_d0.000127613610835.txt",
"cc_n600_k3_d0.000109861142449.txt",
"cc_n700_k3_d9.77987027969e-05.txt",
"cc_n800_k3_d8.98262983103e-05.txt",
"cc_n900_k3_d8.29948594219e-05.txt",
"cc_n1000_k3_d7.68491693592e-05.txt",
"cc_n1100_k3_d7.19527678719e-05.txt",
"cc_n1200_k3_d6.77127865406e-05.txt",
"cc_n1300_k3_d6.40417294615e-05.txt",
"cc_n1400_k3_d5.98609497104e-05.txt",
"cc_n1500_k3_d5.67708081704e-05.txt",
"cc_n1600_k3_d5.44619225457e-05.txt",
"cc_n1700_k3_d5.18041401362e-05.txt",
"cc_n1800_k3_d4.93615853352e-05.txt",
"cc_n1900_k3_d4.76567183783e-05.txt",
"cc_n2000_k3_d4.56033947422e-05.txt",
"cc_n2100_k3_d4.3741099797e-05.txt",
"cc_n2200_k3_d4.17624188026e-05.txt",
"cc_n2300_k3_d4.0385958328e-05.txt",
"cc_n2400_k3_d4.03847906132e-05.txt",
"cc_n2500_k3_d4.03847906132e-05.txt",
"cc_n2600_k3_d3.92668599347e-05.txt",
"cc_n2700_k3_d3.77853568624e-05.txt",
"cc_n2800_k3_d3.61770510155e-05.txt",
"cc_n2900_k3_d3.46813003505e-05.txt",
"cc_n3000_k3_d3.30979158131e-05.txt"
]



def load_hotnet_paths():
    l = []
    for line in hotnet_filenames:
        l.append("../" + network_name + "/out/hotnet2_subnetworks/" + line)
    return l
def load_our_paths(subdir):
    l = []
    our_filenames = ""
    if subdir == "1.0m+1.0c":
        our_filenames = our_filenames_a
    elif subdir == "m*g1cov*g2cov":
        our_filenames = our_filenames_b
    elif subdir == "0.1m+0.9c":
        our_filenames = our_filenames_c
    elif subdir == "0.9m+0.1c":
        our_filenames = our_filenames_d
    elif subdir == "m*g1cov*g2cov+c":
        our_filenames = our_filenames_e
    elif subdir == "m*g1cov*g2cov+c+d":
        our_filenames = our_filenames_f
    for line in our_filenames:
        l.append("../" + network_name + "/out/our_subnetworks/" + subdir + "/" + line)
    return l
