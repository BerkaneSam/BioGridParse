import argparse

parser = argparse.ArgumentParser(description='Process protein interaction databases.\n'
                                             'BioGrid or Intact(filtering not available for IntAct).\n'
                                             'A default search will be done on a BioGrid archetype.')
parser.add_argument('data', nargs='*', metavar='data', type=str, help="The database(s) to be used.")
parser.add_argument('-c', '--codes', nargs='*', help="The different Uniprot codes to be used", required=True)
parser.add_argument('-d', '--document', action='store_true', help="Allows the use of a txt document as argument for"
                                                                  " codes(the file must have uniprot codes separated"
                                                                  " by commas)")
parser.add_argument('-e', '--experience', nargs='*',
                    help="Choose a type of experience to filter with (ex: Two-hybrid, AP-Ms...)")
parser.add_argument('-t', '--throughput', nargs='*',
                    help="Choose a type of throughput to filter with (Low-Throughput,"
                         "High-Throughput, High-Throughput|Low-Throughput)")
parser.add_argument('-i', '--intact', action='store_true', help="Activating parsing on IntAct databases")
parser.add_argument('-a', '--all', action='store_true', help="Activating parsing on all databases available, put "
                                                             "Biogrid database first")
parser.add_argument('-o', '--output', nargs=2, help="naming output files(nodes first, second edges)")
args = parser.parse_args()


def uni_listing(argument):
    """
    Fonction renvoyant les codes uniprot pris en argument sous forme de liste (peut etre ecrit directement ou se
    trouver dans un document .txt au bon format)
    :param argument: code uniprot(-c)
    :return: liste de code uniprot
    """
    print("Registering Uniprot codes to be used...")
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


def parsing(file):
    """
    Fonction permettant de parser la base de donnee BioGrid et de la recuperer sous une forme qui permet l'analyse
    :param file: la base de donnee BIOGRID a parser
    :return: liste de liste contenant les informations de la base de donnee
    """
    print("working on BioGrid...")
    parsedfile = []
    with open(file, "r")as filin:
        for line in filin:
            spl = line.rstrip().split("\t")
            parsedfile.append(spl)
    return parsedfile


def parsing_intact(file):
    """
    Fonction permettant de parser(avec nettoyage des informations inutiles) la base de donnee Intact et de la recuperer
    sous une forme qui permet l'analyse
    :param file: la base de donnee Intact a parser
    :return: liste de liste contenant les informations d'importance de la base de donnee
    """
    print("working on IntAct...")
    parsedfile = []
    with open(file, 'r')as filin:
        for line in filin:
            dataline = []
            if not line.startswith("#"):
                spl = line.rstrip().split("\t")
                '''Gestion d'information perdu dans la base de donnée en ignorant ces lignes (celle avec -)'''
                if spl[1] == '-':
                    continue
                else:
                    int_a = spl[0].split(":")
                    int_b = spl[1].split(":")
                    exp = spl[6].split("(")
                    exp_a = spl[40].split("(")
                    exp_b = spl[41].split("(")
                    pmed = spl[8].split("|")
                    pm = pmed[0].split(":")
                    seq_a = spl[4].split("(")
                    seq_b = spl[5].split("(")
                    seq_aa = seq_a[0].split(":")
                    seq_bb = seq_b[0].split(":")
                    dataline.append(int_a[1])
                    dataline.append(int_b[1])
                    dataline.append(exp[1][0:-1])
                    dataline.append(exp_a[1][0:-1])
                    dataline.append(exp_b[1][0:-1])
                    dataline.append("PUBMED:" + pm[1])
                    dataline.append(seq_aa[1])
                    dataline.append(seq_bb[1])
                    parsedfile.append(dataline)
    return parsedfile


def dico_making(dico, data, i):
    """
    Fonction permettant d'ajouter des elements a un dictionnaire plus specialement ici les informations concernant une
    paire de proteines interagissant ensemble (BioGrid).
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
        '''Gestion de la cle sorted obtenu en verifiant la correspondance par rapport a la base de donnee pour gerer
        l'ordre des autres informations'''
        if key[0] == data[i][23]:
            dico[key]['specie A'].append(data[i][35])
            dico[key]['specie B'].append(data[i][36])
        else:
            dico[key]['specie A'].append(data[i][36])
            dico[key]['specie B'].append(data[i][35])
        dico[key]['getA'].append('NA')
        dico[key]['getB'].append('NA')
        dico[key]['seqA'].append('NA')
        dico[key]['seqB'].append('NA')
    else:
        if key[0] == data[i][23]:
            dc = {'pubmed': [data[i][14]], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                  'specie A': [data[i][35]], 'specie B': [data[i][36]], 'getA': ['NA'], 'getB': ['NA'],
                  'seqA': ['NA'], 'seqB': ['NA']}
            dico[key] = dc
        else:
            dc = {'pubmed': [data[i][14]], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                  'specie A': [data[i][36]], 'specie B': [data[i][35]], 'getA': ['NA'], 'getB': ['NA'],
                  'seqA': ['NA'], 'seqB': ['NA']}
            dico[key] = dc
    return dico


def research(data, codes):
    """
    Fonction permettant de recuperer les experiences/lignes de la base de donnee ou les proteines ont des interactions
    en utilisant les codes Uniprot donnees en arguments et en les comparant a la base de donnee.
    Permet aussi de filtrer par rapport aux types d'experience ou throughput utilise pour mettre en evidence ces
    interactions.
    Recupere aussi les informations a conserver dans un dictionnaire par la suite (BioGrid)
    :param data: base de donnee sous la forme d'une liste obtenu grace a parsing()
    :param codes: liste de code uniprot
    :return: un dictionnaire contenant les informations concernant les interactions entre proteines
    """
    print("parsing of BioGrid data...")
    bd = {}
    for i in range(len(data)):
        if data[i][23] in codes and data[i][26] in codes:
            bd = dico_making(bd, data, i)
    return bd

def predata(data):
    """
    Preparation des donnees de telles sorte que le programme n'est pas un runtime trop élevé
    :param data: donnee a traiter
    :return: liste des données traiter
    """
    print("pretreatment")
    id = []
    prelimid = []
    idkepts = []
    for i in range(len(data)):
        id.append(data[i][5])
    print("part 1 done")
    print(len(id))
    for i in id:
        if i not in prelimid:
            prelimid.append(i)
    print("part 2 done")
    print(len(prelimid))
    for j in range(len(prelimid)):
        #print(id[j])
        print(f"{j} sur {len(prelimid)}")
        if id.count(prelimid[j]) == 1:
            idkepts.append(id[j])
            print(idkepts)
    return idkepts


def explist(data, codes, dico, testid):
    print("experience linking")
    expan = {}
    for i in range(len(data)):
        if data[i][0] in codes and data[i][1] in codes:
            tkey = (data[i][0], data[i][1])
            tkey = tuple(sorted(tkey))
            if tkey in dico and dico[tkey]['pubmed'].count(data[i][5]) == 1 and data[i][5] in testid:
                #print(dico[tkey]['pubmed'])
                for j in range(len(dico[tkey]['pubmed'])):
                    if data[i][5] == dico[tkey]['pubmed'][j]:
                        if dico[tkey]['experiment'][j] not in expan:
                            expan[dico[tkey]['experiment'][j]] = [data[i][2]]
                        elif data[i][2] not in expan[dico[tkey]['experiment'][j]]:
                            expan[dico[tkey]['experiment'][j]].append(data[i][2])
                        else:
                            continue
    return expan

def explisting():
    unicodes = uni_listing(args.codes)
    integretaddata = parsing(args.data[0])
    prelim = research(integretaddata, unicodes)
    del integretaddata
    intdata = parsing_intact(args.data[1])
    testinglist = predata(intdata)
    info = explist(intdata, unicodes, prelim, testinglist)
    for key in info:
        print(f"{key} : {info[key]}")

explisting()





