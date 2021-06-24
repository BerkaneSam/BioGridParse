import argparse

'''
Portion servant à gerer les arguments/parametres du programme
'''
parser = argparse.ArgumentParser(description='getting full list of codes')
parser.add_argument('data', nargs='*', metavar='data', type=str, help="The database(s) to be used.")
parser.add_argument('-b', '--biodata', nargs='*', metavar='data', type=str, help="biogrid list codes", required=True)
args = parser.parse_args()


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


def getnames(data):
    """
    Fonction permettant de recuperer tous les codes proteiques present dans la base de donnee apres traitement sous
    forme d'une liste tout en evitant la presence de doublon
    :param data: base de donnee apres traitement par parsing_intact()
    :return: liste de codes proteiques
    """
    print("working on IntAct codes...")
    nameslist = []
    for i in range(len(data)):
        if data[i][0] not in nameslist:
            nameslist.append(data[i][0])
        if data[i][1] not in nameslist:
            nameslist.append(data[i][1])
    return nameslist


def namefile(data, mainfile):
    """
    Fonction permettant de creer un fichier txt contenant tous les codes IntAct et Biogrid en liant les codes recuperer
    à l'aide du programme parsebio.py a ceux recuperer par getnames()
    :param data: liste de codes proteiques
    :param mainfile: fichier contenant les codes proteiques BioGrid
    :return: un fichier txt avec tous les codes proteiques present dans Biogrid et IntAct
    """
    print("working on BioGrid and Intact codes...")
    biolist = []
    with open(mainfile, 'r') as filin:
        for line in filin:
            spl = line.rstrip().split(",")
            for i in spl:
                biolist.append(i)
    for i in data:
        if i not in biolist:
            biolist.append(i)
    print("making a file...")
    with open("intbiolistcodes.txt", 'w') as fileout:
        for i in biolist:
            i = str(i)
            fileout.write(f"{i},")


def fulllaunch():
    """
    Programme permettant d'obtenir le fichier txt en utilisant toutes les fonctions pour BioGrid et IntAct
    :return: un fichier txt
    """
    print("program launched")
    intactdata = parsing_intact(args.data[0])
    codelist = getnames(intactdata)
    namefile(codelist, args.biodata[0])
    print("processing finished")


fulllaunch()
