import argparse

'''
Portion servant à gerer les arguments/parametres du programme
'''
parser = argparse.ArgumentParser(description='getting full list of codes')
parser.add_argument('data', nargs='*', metavar='data', type=str, help="The database(s) to be used.")
args = parser.parse_args()

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

def getcodes(data):
    print("working on codes...")
    codelist = []
    for i in range(len(data)):
        if data[i][23] not in codelist:
            codelist.append(data[i][23])
        if data[i][26] not in codelist:
            codelist.append(data[i][26])
    return codelist

def filemaking(data):
    print("making file...")
    with open("biolistcodes.txt",'w') as fileout:
        for i in data:
            i = str(i)
            fileout.write(f"{i},")

def launch():
    biodata = parsing(args.data[0])
    codelist = getcodes(biodata)
    filemaking(codelist)

launch()