import argparse

'''
Portion servant à gerer les arguments du programme
'''
parser = argparse.ArgumentParser(description='Process BioGrid databases.')
parser.add_argument('data', metavar='data', type=str, help="The BIOGGRID database to be used.")
parser.add_argument('-c', '--codes', nargs=2, help="The different Uniprot codes to be used", required=True)
parser.add_argument('-e', '--experience', help="Choose a type of experience to filter with (ex: Two-hybrid, AP-Ms...)")
parser.add_argument('-t', '--throughput',
                    help="Choose a type of throughput to filter with (Low-Throughput,"
                         "High-Throughput, High-Throughput|Low-Throughput)")
args = parser.parse_args()


def switch_exp(argument):
    """
    Fonction permettant de traduire les arguments lie aux experiences vers leurs versions dans la base de donnee
    argument : l'argument d'experience(-e)
    :return: un champ txt avec l'argument traduit ou message par défaut si absence d'argument
    """
    switcher = {
        "AP-MS": 'Affinity Capture-MS',
        "Two-hybrid": 'Two-hybrid',
        "RC": 'Reconstituted Complex',
        "ACW": 'Affinity Capture-Western',
        "ACL": 'Affinity Capture-Luminescence',
        "BA": 'Biochemical Activity',
        "CCS": 'Co-crystal Structure',
        "FW": 'Far Western',
        "FRET": 'FRET',
        "Protein-peptide": 'Protein-peptide',
        "Co-localization": 'Co-localization',
        "AC-RNA": 'Affinity Capture-RNA',
        "Protein-RNA": 'Protein-RNA',
        "PCA": 'PCA',
        "Co-purification": 'Co-purification',
        "Co-fractionation": 'Co-fractionation',
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
        "Low-Throughput": 'Low Throughput',
        "High-Throughput": 'High Throughput',
        "High-Throughput|Low-Throughput": 'High Throughput|Low Throughput'
    }
    return switcher.get(argument, "No filter")


def parsing(file):
    """
    Fonction permettant de parser la base de donnees et de la recuperer sous une forme qui permet l'analyse
    :param file: la base de donnee BIOGRID a parser
    :return: liste de liste contenant les informations de la base de donnee
    """
    parsedfile = []
    with open(file, "r")as filin:
        for line in filin:
            line.rstrip()
            spl = line.split("\t")
            parsedfile.append(spl)
    return parsedfile


def research(data):
    """
    Fonction permettant de recuperer les experiences/lignes de la base de donnee ou les proteines ont des interactions
    en utilisant les codes Uniprot donnees en arguments et en les comparant a la base de donnee.
    Permet aussi de filtrer par rapport aux types d'experience utilise pour mettre en evidence ces interactions.
    :param data: base de donnee sous la forme d'une liste obtenu grace a parsing().
    :return: une liste contenant les numero de ligne de la base de donnee ou les conditions sont rempli
    """
    dataline = []
    for i in range(len(data)):
        if (data[i][23] != args.codes[0] or data[i][26] != args.codes[1]) and (
                data[i][23] != args.codes[1] or data[i][26] != args.codes[0]):
            continue
        if switch_exp(args.experience) != "Invalid experiment":
            if data[i][11] == switch_exp(args.experience):
                dataline.append(i + 1)
        else:
            dataline.append(i + 1)
    return dataline


def research_through(vector, data):
    """
    Une fonction permettant de tester sous quelle throughput l'experience a ete faite si mis en argument
    :param vector: liste contenant les resultat preliminaire obtenu a l'aide de research()
    :param data: la base de donnee sous forme de liste obtenu a l'aide de parsing()
    :return: une liste contenant les numero de ligne de la base de donnee ou les conditions sont rempli
    """
    dataline = []
    if switch_throu(args.throughput) != "No filter":
        for i in vector:
            if data[i - 1][17] == switch_throu(args.throughput):
                dataline.append(i)
        return dataline
    else:
        return vector


def biogrind():
    """
    Le programme utilisant toutes les fonctions permettant l'analyse de la base de donnees souhaiter
    :return: impression des resultat
    """
    integretaddata = parsing(args.data)
    preliminary = research(integretaddata)
    occurences = research_through(preliminary, integretaddata)
    if not occurences:
        print("Pas d'experience rassemblant ces protéines (possiblement liées aux conditions)")
    else:
        print(f"Lignes ou l'on retrouve ces deux codes Uniprot : {occurences}")


biogrind()
