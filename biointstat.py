import argparse
import matplotlib.pyplot as plt

'''
Portion servant à gerer les arguments/parametres du programme
'''
parser = argparse.ArgumentParser(description='Process a protein network to find the shortest path between two proteins')
parser.add_argument('network', metavar='network', type=str,
                    help="The network to be used(edges).")
parser.add_argument('-s', '--shift', nargs=1, help="choosing an interval for the barplot", required=True)
args = parser.parse_args()


def parse_edges(edges):
    """
    Fonction permettant le traitement des données du fichier JSON contenant les edges afin de récuperer toutes les
    interactions(edges) et de renvoyer une liste de ces edges
    :param edges: fichier JSON des edges
    :return: liste des edges
    """
    print("parsing edges...")
    parsed_edges = []
    with open(edges, 'r', encoding='utf8') as filin:
        for line in filin:
            spl = line.split("'")
            pre_edges = ([spl[1], spl[3]])
            pre_edges = tuple(sorted(pre_edges))
            if pre_edges[0] != pre_edges[1]:
                parsed_edges.append(pre_edges)
    "Ce passage par un dictionnaire permet d'eviter toute redondance"
    parsed_edges = list(dict.fromkeys(parsed_edges))
    return parsed_edges


def inter_prot(edges):
    """
    Fonction permettant de trouver le nombre d'interaction de chaque noeud (ici, protéine) et de garder ce comptage dans
    un dictionnaire
    :param edges: les edges traiter précedemment à l'aide de la fonction parse_edges()
    :return: un dictionnaire contenant le comptage pour chaque protéine
    """
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
    """
    Cette fonction permet de recuperer un dictionnaire ayant pour cle des intervalles de valeurs representant le nombre
    d'interaction (edges) d'une proteine et ces cle ont pour valeur un comptage du nombre de proteine ayant des valeurs
    d'interactions presente dans cette intervalle
    :param pdata: un dictionnaire de protéine avec leurs nombres d'interactions obtenu avec inter_prot()
    :param shift: l'espace des intervalles
    :param maxlen: la plus haute valeur d'interactions d'une protéine présente dans pdata obtenu a l'aide de highest()
    :return: un dictionnaire avec pour cle des intervalles de nombres d'interactions et pour valeur le nombre de
    proteine avec un nombre d'interaction compris dans cette intervalle
    """
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
    """
    Fonction permettant d'obtenir la valeur d'interactions maximale d'une proteine presente dans le dictionnaire
    :param dico: un dictionnaire de protéine avec leurs nombres d'interactions obtenu avec inter_prot()
    :return: la valeur d'interaction maximale
    """
    higher = 0
    for key in dico:
        tmp = dico[key]
        if tmp > higher:
            higher = tmp
    return higher


def stat_inter():
    """
    Le programme utilisant toutes les fonctions et permettant à l'aide de leurs resultat d'afficher un histogramme
    :return: un histogramme
    """
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
