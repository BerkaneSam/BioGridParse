import argparse

'''
Portion servant à gerer les arguments du programme
'''
parser = argparse.ArgumentParser(description='Process BioGrid databases.')
parser.add_argument('data', metavar='data', type=str, help="The BIOGGRID database to be used.")
parser.add_argument('-c', '--codes', nargs='*', help="The different Uniprot codes to be used", required=True)
parser.add_argument('-d', '--document', action='store_true', help="Allows the use of a txt document as argument for"
                                                                  " codes(the file must have uniprot codes separated"
                                                                  " by commas)")
parser.add_argument('-e', '--experience', nargs='*',
                    help="Choose a type of experience to filter with (ex: Two-hybrid, AP-Ms...)")
parser.add_argument('-t', '--throughput', nargs='*',
                    help="Choose a type of throughput to filter with (Low-Throughput,"
                         "High-Throughput, High-Throughput|Low-Throughput)")
args = parser.parse_args()


def uni_listing(argument):
    """
    Fonction renvoyant les codes uniprot pris en argument sous forme de liste
    :param argument: code uniprot(-c)
    :return: liste de code uniprot
    """
    unicodes = []
    if args.document:
        with open(argument[0], 'r')as filin:
            for line in filin:
                rline = line.rstrip()
                spl = rline.split(",")
                for i in spl:
                    unicodes.append(i)
    else:
        for i in argument:
            if i not in unicodes:
                unicodes.append(i)
    return unicodes


def switch_exp(argument):
    """
    Fonction permettant de traduire les arguments lie aux experiences vers leurs versions dans la base de donnee
    argument : l'argument d'experience(-e)
    :return: un champ txt avec l'argument traduit ou message par défaut si absence d'argument
    """
    switcher = {
        "AP-MS": 'Affinity Capture-MS',
        "T-H": 'Two-hybrid',
        "RC": 'Reconstituted Complex',
        "ACW": 'Affinity Capture-Western',
        "ACL": 'Affinity Capture-Luminescence',
        "BA": 'Biochemical Activity',
        "CCS": 'Co-crystal Structure',
        "FW": 'Far Western',
        "FRET": 'FRET',
        "PP": 'Protein-peptide',
        "Co-l": 'Co-localization',
        "AC-RNA": 'Affinity Capture-RNA',
        "P-RNA": 'Protein-RNA',
        "PCA": 'PCA',
        "Co-p": 'Co-purification',
        "Co-f": 'Co-fractionation',
        "PL-MS": 'Proximity Label-MS'
    }
    return switcher.get(argument, "Invalid experiment")


def switch_throu(argument):
    """
    Fonction permettant de traduire les arguments lie aux throughput vers leurs versions dans la base de donnee
    :param argument: l'argument de throughput (-t)
    :return: un champ txt avec l'argument traduit ou message par défaut si absence d'argument
    """
    switcher = {
        "LT": 'Low Throughput',
        "HT": 'High Throughput',
        "HLT": 'High Throughput|Low Throughput'
    }
    return switcher.get(argument, "No filter")


def filter_listing(argument, ty):
    """
    Fonction permettant de renvoyer une liste d'element avec lesquelle filtrer la recherche
    :param argument: argument de filtrage
    :param ty: type de filtrage (e: experiment, t: throughput
    :return: return un liste de terme qui serviront de filtre
    """
    filt = []
    if argument:
        if ty == 'e':
            for i in argument:
                tr = switch_exp(i)
                if tr not in filt:
                    filt.append(tr)
        if ty == 't':
            for i in argument:
                tr = switch_throu(i)
                if tr not in filt:
                    filt.append(tr)
    return filt


def filter_error(argument1, argument2):
    """
    Gestion de mauvaises informations données en filtrage
    :param argument1: filtrages d'experience
    :param argument2: filtrage de throughput
    :return: valeurs booleene
    """
    n = 0
    m = 0
    if argument1:
        for i in range(len(argument1)):
            if switch_exp(argument1[i]) == "Invalid experiment":
                n = 1
    if argument2:
        for i in range(len(argument2)):
            if switch_throu(argument2[i]) == "No filter":
                m = 1
    return n, m


def parsing(file):
    """
    Fonction permettant de parser la base de donnees et de la recuperer sous une forme qui permet l'analyse
    :param file: la base de donnee BIOGRID a parser
    :return: liste de liste contenant les informations de la base de donnee
    """
    print("working on data...")
    parsedfile = []
    with open(file, "r")as filin:
        for line in filin:
            line.rstrip()
            spl = line.split("\t")
            parsedfile.append(spl)
    return parsedfile


def dico_making(dico, data, i):
    """
    Fonction permettant d'ajouter des elements a un dictionnaire plus specialement ici les informations concernant une
    paire de proteines interagissant ensemble.
    :param dico: dictionnaire sur lequel la fonction travaille
    :param data: le jeu de donnees ou les informations sont recupere
    :param i: parametre de position
    :return: un dictionnaire avec des ajouts
    """
    key = (data[i][23], data[i][26])
    key = tuple(sorted(key))
    if key in dico:
        dico[key]['pubmed'].append(data[i][14])
        dico[key]['experiment'].append(data[i][11])
        dico[key]['throughput'].append(data[i][17])
    else:
        dc = {'pubmed': [data[i][14]], 'experiment': [data[i][11]], 'throughput': [data[i][17]]}
        dico[key] = dc
    return dico


def research(data, codes):
    """
    Fonction permettant de recuperer les experiences/lignes de la base de donnee ou les proteines ont des interactions
    en utilisant les codes Uniprot donnees en arguments et en les comparant a la base de donnee.
    Permet aussi de filtrer par rapport aux types d'experience ou throughput utilise pour mettre en evidence ces
    interactions.
    Recupere aussi les informations a conserver dans un dictionnaire par la suite
    :param data: base de donnee sous la forme d'une liste obtenu grace a parsing().
    :param codes: liste de code uniprot
    :return: un dictionnaire contenant les informations concernant les interactions entre deux proteines
    """
    print("parsing of the data...")
    bd = {}
    exp = filter_listing(args.experience, 'e')
    throu = filter_listing(args.throughput, 't')
    for i in range(len(data)):
        if data[i][23] in codes and data[i][26] in codes:
            if (exp and "Invalid experiment" not in exp) and (throu and "No filter" not in throu):
                if data[i][11] in exp and data[i][17] in throu:
                    bd = dico_making(bd, data, i)
            elif exp and "Invalid experiment" not in exp:
                if data[i][11] in exp:
                    bd = dico_making(bd, data, i)
            elif throu and "No filter" not in throu:
                if data[i][17] in throu:
                    bd = dico_making(bd, data, i)
            else:
                bd = dico_making(bd, data, i)
    return bd


def node_making_json(codes):
    """
    Fonction permettant de creer les nodes au format Json, en creant un node pour chaque code uniprot en parametre
    :param codes: les codes uniprot donnees en parametres recupere sous forme de liste a l'aide de la fonction
    uni_listing
    :return: imprime les nodes (dans un fichier a l'aide de >)
    """
    print("Making nodes.txt...")
    with open("nodes.txt", 'w')as filout:
        for i in codes:
            i = str(i)
            filout.write(f'{{ data: {{ id: \'{i}\', size: 200, name: \'{i}\' }}, classes: [] }},\n')


def edge_making_json(dico):
    """
    Fonction permettant de creer les edges sous format Json, en creant un edge pour chaque paire de proteine en
    interaction.
    Les informations concernant les paires de proteines, ainsi que les paires de proteines sont recupere depuis
    dico.
    :param dico: dictionnaire contenant les informations concernant les proteines en interactions
    :return: imprime les edges (dans un fichier a l'aide de >)
    """
    print("Making edge.txt...")
    with open("edge.txt", 'w')as filout:
        for key in dico:
            classes = []
            for i in range(len(dico[key]['pubmed'])):
                classes.append(dico[key]['pubmed'][i])
                classes.append(dico[key]['experiment'][i])
                classes.append(dico[key]['throughput'][i])
            filout.write(
                f"{{ data: {{ source: \'{key[0]}\', target: \'{key[1]}\', strength: 5 }}, classes: {classes} }},\n")


def biogrind():
    """
    Le programme utilisant toutes les fonctions permettant l'analyse de la base de donnees souhaiter
    :return: impression des resultat
    """
    print("Program launched")
    exptest, throutest = filter_error(args.experience, args.throughput)
    if exptest == 1:
        print("Une erreur est survenu, lie aux filtre des experiences verifier que vous respectez les bonnes"
              " conditions!")
        exit()
    if throutest == 1:
        print("Une erreur est survenu, lie aux filtre des throughput verifier que vous respectez les bonnes"
              " conditions!")
        exit()
    integretaddata = parsing(args.data)
    unicodes = uni_listing(args.codes)
    occurences = research(integretaddata, unicodes)
    if not occurences:
        print("Pas de correspondance avec la recherche trouve")
        exit()
    node_making_json(unicodes)
    edge_making_json(occurences)


biogrind()
