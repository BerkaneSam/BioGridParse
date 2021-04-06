def parsing(file):
    parsedfile = []
    with open(file, "r")as filin:
        for line in filin:
            line.rstrip()
            spl = line.split("\t")
            parsedfile.append(spl)
    return parsedfile


def research(data):
    dataline = []
    print("Entrer le premier code Uniprot : ")
    first = input()
    print("Entrer le second code Uniprot : ")
    second = input()
    for i in range(len(data)):
        if (data[i][23] != first or data[i][26] != second) and (data[i][23] != second or data[i][26] != first):
            continue
        dataline.append(i + 1)
    return dataline


def biogrind():
    # print("Enter the BioGrid database to search on : ")
    # database = input()
    database = "BIOGRID-MV-Physical-4.3.196.tab3.txt"
    integretaddata = parsing(database)
    occurences = research(integretaddata)
    if not occurences:
        print("Pas d'experience rassemblant ces protéines")
    else:
        print(f"Lignes ou l'on retrouve ces deux codes Uniprot : {occurences}")


biogrind()
