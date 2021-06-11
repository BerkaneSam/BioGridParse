import argparse

'''
Portion servant à gerer les arguments/parametres du programme
'''
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

'''
Fonction utilise par le programme
'''


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
            pm = None
            doi = None
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
                    for i in range(len(pmed)):
                        if pmed[i].startswith("pub"):
                            pm = pmed[i].split(":")
                        if pmed[i].startswith("doi"):
                            doi = pmed[i].split(":")
                    # pm = pmed[0].split(":")
                    seq_a = spl[4].split("(")
                    seq_b = spl[5].split("(")
                    seq_aa = seq_a[0].split(":")
                    seq_bb = seq_b[0].split(":")
                    dataline.append(int_a[1])
                    dataline.append(int_b[1])
                    dataline.append(exp[1][0:-1])
                    dataline.append(exp_a[1][0:-1])
                    dataline.append(exp_b[1][0:-1])
                    if pm is not None:
                        dataline.append("PUBMED:" + pm[1])
                    else:
                        dataline.append("NA")
                    dataline.append(seq_aa[1])
                    dataline.append(seq_bb[1])
                    if doi is not None:
                        dataline.append("DOI:" + doi[1])
                    else:
                        dataline.append("NA")
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
        if data[i][14].startswith("PUB"):
            dico[key]['pubmed'].append(data[i][14])
            dico[key]['DOI'].append('NA')
        else:
            dico[key]['pubmed'].append('NA')
            dico[key]['DOI'].append(data[i][14])
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
            if data[i][14].startswith("PUB"):
                dc = {'pubmed': [data[i][14]], 'DOI': ['NA'], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                      'specie A': [data[i][35]], 'specie B': [data[i][36]], 'getA': ['NA'], 'getB': ['NA'],
                      'seqA': ['NA'], 'seqB': ['NA']}
                dico[key] = dc
            else:
                dc = {'pubmed': ['NA'], 'DOI': [data[i][14]], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                      'specie A': [data[i][35]], 'specie B': [data[i][36]], 'getA': ['NA'], 'getB': ['NA'],
                      'seqA': ['NA'], 'seqB': ['NA']}
                dico[key] = dc
        else:
            if data[i][14].startswith("PUB"):
                dc = {'pubmed': [data[i][14]], 'DOI': ['NA'], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                      'specie A': [data[i][36]], 'specie B': [data[i][35]], 'getA': ['NA'], 'getB': ['NA'],
                      'seqA': ['NA'], 'seqB': ['NA']}
                dico[key] = dc
            else:
                dc = {'pubmed': ['NA'], 'DOI': [data[i][14]], 'experiment': [data[i][11]], 'throughput': [data[i][17]],
                      'specie A': [data[i][36]], 'specie B': [data[i][35]], 'getA': ['NA'], 'getB': ['NA'],
                      'seqA': ['NA'], 'seqB': ['NA']}
                dico[key] = dc
    return dico


def dico_making_intact(dico, data, i):
    """
    Fonction permettant d'ajouter des elements a un dictionnaire plus specialement ici les informations concernant une
    paire de proteines interagissant ensemble (Intact).
    :param dico: dictionnaire sur lequel la fonction travaille
    :param data: le jeu de donnees ou les informations sont recupere
    :param i: parametre de position
    :return: un dictionnaire avec des ajouts
    """
    key = (data[i][0], data[i][1])
    key = tuple(sorted(key))
    if key in dico:
        dico[key]['pubmed'].append(data[i][5])
        dico[key]['DOI'].append(data[i][8])
        dico[key]['experiment'].append(data[i][2])
        dico[key]['throughput'].append('NA')
        dico[key]['specie A'].append('NA')
        dico[key]['specie B'].append('NA')
        if key[0] == data[i][0]:
            dico[key]['getA'].append(data[i][3])
            dico[key]['getB'].append(data[i][4])
            dico[key]['seqA'].append(data[i][6])
            dico[key]['seqB'].append(data[i][7])
        else:
            dico[key]['getA'].append(data[i][4])
            dico[key]['getB'].append(data[i][3])
            dico[key]['seqA'].append(data[i][7])
            dico[key]['seqB'].append(data[i][6])
    else:
        if key[0] == data[i][0]:
            dc = {'pubmed': [data[i][5]], 'DOI': [data[i][8]], 'experiment': [data[i][2]], 'throughput': ['NA'],
                  'specie A': ['NA'],
                  'specie B': ['NA'], 'getA': [data[i][3]], 'getB': [data[i][4]], 'seqA': [data[i][6]],
                  'seqB': [data[i][7]]}
            dico[key] = dc
        else:
            dc = {'pubmed': [data[i][5]], 'DOI': [data[i][8]], 'experiment': [data[i][2]], 'throughput': ['NA'],
                  'specie A': ['NA'],
                  'specie B': ['NA'], 'getA': [data[i][4]], 'getB': [data[i][3]], 'seqA': [data[i][7]],
                  'seqB': [data[i][6]]}
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


def research_intact(data, codes, bd=None):
    """
    Fonction permettant de recuperer les experiences/lignes de la base de donnee ou les proteines ont des interactions
    en utilisant les codes Uniprot donnees en arguments et en les comparant a la base de donnee.
    Recupere aussi les informations a conserver dans un dictionnaire par la suite (Intact)
    :param data: base de donnee sous la forme d'une liste obtenu grace a parsing_intact()
    :param codes: liste de code uniprot
    :param bd: dico a utiliser si un dico deja rempli disponible
    :return: un dictionnaire contenant les informations concernant les interactions entre proteines
    """
    print("parsing of IntAct data...")
    if bd is None:
        bd = {}
    for i in range(len(data)):
        if data[i][0] in codes and data[i][1] in codes:
            bd = dico_making_intact(bd, data, i)
    return bd


def node_making_json(codes, output="nodes.txt"):
    """
    Fonction permettant de creer les nodes au format Json, en creant un node pour chaque code uniprot en parametre
    :param codes: les codes uniprot donnees en parametres recupere sous forme de liste a l'aide de la fonction
    uni_listing
    :param output: nom du fichier de sortie
    :return: imprime les nodes (dans un fichier a l'aide de >)
    """
    print("Making nodes.txt...")
    with open(output, 'w')as filout:
        for i in codes:
            i = str(i)
            filout.write(f'{{ data: {{ id: \'{i}\', size: 200, name: \'{i}\' }}, classes: [] }},\n')


def edge_making_json(dico, output="edges.txt"):
    """
    Fonction permettant de creer les edges sous format Json, en creant un edge pour chaque paire de proteine en
    interaction.
    Les informations concernant les paires de proteines, ainsi que les paires de proteines sont recupere depuis
    dico.
    :param dico: dictionnaire contenant les informations concernant les proteines en interactions
    :param output: nom du fichier de sortie
    :return: imprime les edges (dans un fichier a l'aide de >)
    """
    print("Making edge.txt...")
    with open(output, 'w')as filout:
        for key in dico:
            classes = []
            for i in range(len(dico[key]['pubmed'])):
                classes.append(dico[key]['pubmed'][i])
                classes.append(dico[key]['DOI'][i])
                classes.append(dico[key]['experiment'][i])
                classes.append(dico[key]['throughput'][i])
                classes.append(dico[key]['specie A'][i])
                classes.append(dico[key]['specie B'][i])
                classes.append(dico[key]['getA'][i])
                classes.append(dico[key]['getB'][i])
                classes.append(dico[key]['seqA'][i])
                classes.append(dico[key]['seqB'][i])
            filout.write(
                f"{{ data: {{ source: \'{key[0]}\', target: \'{key[1]}\', strength: 5 }}, classes: {classes} }},\n")


def coverage(data, codes, dico, size):
    print("coverage calculated...")
    oncount = []
    for i in range(len(data)):
        if data[i][0] in codes and data[i][1] in codes:
            tkey = (data[i][0], data[i][1])
            tkey = tuple(sorted(tkey))
            if tkey in dico and data[i][5] in dico[tkey]['pubmed']:
                oncount.append(tkey)
    print(f"percentage coverage of BioGrid {(len(oncount) * 100) / size}%")
    print(f"percentage coverage of IntAct {(len(oncount) * 100) / len(data)}%")
    print(f"number of interactions found one both databases : {len(oncount)}")


'''
Le programme
'''


def biogrind():
    """
    Le programme utilisant toutes les fonctions permettant l'analyse de la base de donnees souhaite
    :return: impression des resultat dans deux fichiers nodes.txt et edge.txt au format JSON pour cytoscape
    """
    print("Program launched")
    '''Gestion de deux options incompatibles'''
    if args.intact and args.all:
        print("error : parameters -i and -a can't be used at the same time")
        exit()
    if args.intact:
        '''Lancement sur une base de donnee IntAct'''
        print("no filters available for IntAct")
        intdata = parsing_intact(args.data[0])
        unicodes = uni_listing(args.codes)
        occurences = research_intact(intdata, unicodes)
    elif args.all:
        '''Lancement sur base de donnee BioGrid et IntAct, parcourt BioGrid en premier'''
        print("filters only used on BioGrid database")
        unicodes = uni_listing(args.codes)
        integretaddata = parsing(args.data[0])
        prelim = research(integretaddata, unicodes)
        # bio_size = len(integretaddata) - 1
        del integretaddata
        intdata = parsing_intact(args.data[1])
        # coverage(intdata, unicodes, prelim, bio_size)
        occurences = research_intact(intdata, unicodes, prelim)
        del intdata
    else:
        '''Lancement par defaut sur une base de donnee BioGrid'''
        exptest, throutest = filter_error(args.experience, args.throughput)
        '''Gestion des erreurs lie au filtrages'''
        if exptest == 1:
            print("Une erreur est survenu, lie aux filtre des experiences verifier que vous respectez les bonnes"
                  " conditions!")
            exit()
        if throutest == 1:
            print("Une erreur est survenu, lie aux filtre des throughput verifier que vous respectez les bonnes"
                  " conditions!")
            exit()
        integretaddata = parsing(args.data[0])
        unicodes = uni_listing(args.codes)
        occurences = research(integretaddata, unicodes)
    '''Gestion d'absence de resultat'''
    if not occurences:
        print("Pas de correspondance avec la recherche trouve")
        exit()
    '''Creation des fichiers JSON contenant les nodes et edges'''
    if args.output:
        node_making_json(unicodes, args.output[0])
        edge_making_json(occurences, args.output[1])
    else:
        node_making_json(unicodes)
        edge_making_json(occurences)
    print("Program done")


'''
Lancement du programme
'''
biogrind()
