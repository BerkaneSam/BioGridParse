import argparse

parser = argparse.ArgumentParser(description='Process a protein network to find the shortest path between two proteins')
parser.add_argument('network', nargs=2, metavar='network', type=str,
                    help="The network to be used(first nodes then edges).")
parser.add_argument('-c', '--codes', nargs=2, help="The different Uniprot codes to be used", required=True)
args = parser.parse_args()


def parse_nodes(nodes):
    parsed_nodes = {}
    with open(nodes, 'r') as filin:
        for line in filin:
            spl = line.rstrip().split("'")
            print(spl)
            parsed_nodes[spl[1]] = False
    return parsed_nodes


def parse_edges(edges):
    parsed_edges = []
    with open(edges, 'r') as filin:
        for line in filin:
            pre_edges = []
            spl = line.split("'")
            pre_edges.extend([spl[1], spl[3]])
            if pre_edges[0] != pre_edges[1]:
                parsed_edges.append(pre_edges)
    return parsed_edges


def uni_listing(argument):
    unicodes = []
    for i in argument:
        if i not in unicodes:
            unicodes.append(i)
    return unicodes


def check_nodes(nodes, codes):
    for i in codes:
        if i not in nodes.keys():
            return 1
    return 0


def shortest_path(inter, main, target, track, marked):
    queue = [main]
    marked[main] = True
    while len(queue) != 0:
        x = queue[0]
        queue.pop(0)
        for edge in inter:
            if x in edge:
                for node in edge:
                    if node != x and marked[node] is not True:
                        queue.append(node)
                        track[node] = x
                        marked[x] = True
                        if node == target:
                            return True
    return False


def find_path(inter, main, target, marked):
    tracker = {}
    link = []
    if not shortest_path(inter, main, target, tracker, marked):
        print("The two proteins aren't linked")
    else:
        link.append(target)
        pathing = target
        while tracker[pathing] != main:
            link.append(tracker[pathing])
            pathing = tracker[pathing]
        link.append(main)
    return link


def shortpath_main():
    nodes = parse_nodes(args.network[0])
    unicodes = uni_listing(args.codes)
    check = check_nodes(nodes, unicodes)
    if check == 1:
        print("At least one of the proteins is not present")
        exit()
    edges = parse_edges(args.network[1])
    shorter_path = find_path(edges, args.codes[0], args.codes[1], nodes)
    print(shorter_path)


shortpath_main()
