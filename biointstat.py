import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Process a protein network to find the shortest path between two proteins')
parser.add_argument('network', metavar='network', type=str,
                    help="The network to be used(edges).")
parser.add_argument('-s', '--shift', nargs=1, help="choosing an interval for the barplot", required=True)
args = parser.parse_args()


def parse_edges(edges):
    print("parsing edges...")
    parsed_edges = []
    with open(edges, 'r', encoding='utf8') as filin:
        for line in filin:
            spl = line.split("'")
            pre_edges = ([spl[1], spl[3]])
            pre_edges = tuple(sorted(pre_edges))
            if pre_edges[0] != pre_edges[1]:
                parsed_edges.append(pre_edges)
    parsed_edges = list(dict.fromkeys(parsed_edges))
    return parsed_edges


def inter_prot(edges):
    print("counting interactions...")
    inter = {}
    for edge in edges:
        if edge[0] in inter:
            inter[edge[0]] += 1
        else:
            inter[edge[0]] = 1
        if edge[1] in inter:
            inter[edge[1]] += 1
        else:
            inter[edge[1]] = 1
    return inter


def new_data(pdata, shift, maxlen):
    print("data treatment...")
    shift = int(shift)
    base = 0
    top = shift
    intervals = []
    data = {}
    print(" filling empty dict...")
    while base < maxlen:
        cur_key = (base, top)
        intervals.append(cur_key)
        base += shift
        top += shift
    print(" filling dict values...")
    for interval in intervals:
        for key in pdata:
            if interval[0] < pdata[key] <= interval[1]:
                if interval[1] in data:
                    data[interval[1]] += 1
                else:
                    data[interval[1]] = 1
    return data


def highest(dico):
    higher = 0
    for key in dico:
        tmp = dico[key]
        if tmp > higher:
            higher = tmp
    return higher


def stat_inter():
    print("launching program")
    data = parse_edges(args.network)
    result = inter_prot(data)
    for key in result:
        print(result[key])
    print("done counting")
    maxlen = highest(result)
    fdata = new_data(result, args.shift[0], maxlen)
    print(maxlen)
    plt.bar(list(fdata.keys()), fdata.values(), color='r')
    plt.title("Répartition des protéines selon leurs nombres d'edges")
    plt.xlabel("Intervalle représentant le nombre d'edge")
    plt.ylabel("Nombre de protéines par intervalle d'edge")
    plt.show()
    print("program done")

stat_inter()
