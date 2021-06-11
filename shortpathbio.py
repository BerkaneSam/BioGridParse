import argparse

parser = argparse.ArgumentParser(description='Process a protein network to find the shortest path between two proteins')
parser.add_argument('network', nargs=2, metavar='network', type=str,
                    help="The network to be used(first nodes then edges).")
parser.add_argument('-c', '--codes', nargs=2, help="The different Uniprot codes to be used", required=True)
parser.add_argument('-n', '--neighbors', action='store_true', help="finding neighboring interactions of proteins"
                                                                   "in the shortest path interactions")
parser.add_argument('-o', '--output', nargs=2, help="First name nodes output then edges")
args = parser.parse_args()


def parse_nodes(nodes):
    print("parsing nodes...")
    parsed_nodes = {}
    with open(nodes, 'r') as filin:
        for line in filin:
            spl = line.rstrip().split("'")
            parsed_nodes[spl[1]] = False
    return parsed_nodes


def parse_edges(edges):
    print("parsing edges...")
    parsed_edges = []
    edges_depository = []
    with open(edges, 'r') as filin:
        for line in filin:
            edges_depository.append(line)
            pre_edges = []
            spl = line.split("'")
            pre_edges.extend([spl[1], spl[3]])
            if pre_edges[0] != pre_edges[1]:
                parsed_edges.append(pre_edges)
    edges_depository = list(dict.fromkeys(edges_depository))
    return parsed_edges, edges_depository


def uni_listing(argument):
    print("code listing...")
    unicodes = []
    for i in argument:
        if i not in unicodes:
            unicodes.append(i)
    return unicodes


def check_nodes(nodes, codes):
    print("checking nodes...")
    for i in codes:
        if i not in nodes.keys():
            return 1
    return 0


def shortest_path(inter, main, target, track, marked):
    print("BFS ongoing...")
    queue = [main]
    marked[main] = True
    dist = {main: 0}
    tmp = None
    while len(queue) != 0:
        x = queue[0]
        queue.pop(0)
        for edge in inter:
            if x in edge:
                for node in edge:
                    if marked[node] is not True:
                        dist[node] = dist[x] + 1
                        queue.append(node)
                        marked[node] = True
                        if node == target:
                            tmp = dist[node]
                    if node != x:
                        if tmp is not None and dist[node] > tmp:
                            return True, tmp
                        if dist[x] < dist[node]:
                            if node in track.keys() and x not in track[node]:
                                track[node].append(x)
                            else:
                                track[node] = [x]
    return False


def pathfinding(path, main, target, limit):
    mpath = []
    if main == target:
        return [[main]]
    for node in path[target]:
        cur_node = pathfinding(path, main, node, limit - 1)
        for lnodes in cur_node:
            if len(lnodes)+1 < limit:
                lnodes.append(target)
                mpath.append(lnodes)
    return mpath


def find_path(inter, main, target, marked):
    print("launch of shortest path...")
    tracker = {}
    link = []
    result, dist = shortest_path(inter, main, target, tracker, marked)
    print(dist)
    if not result:
        print("The two proteins aren't linked")
    else:
        link = pathfinding(tracker, main, target, dist + 2)
        print(link)
    return link


def neighbor_finder(path, edges):
    print("looking for neighbors...")
    neighbors = {}
    mem = []
    for code in path:
        print(f"{path.index(code)} sur {len(path)}")
        for edge in edges:
            if code in edge:
                for node in edge:
                    if node != code and node not in path:
                        if code in neighbors:
                            neighbors[code].append(node)
                        else:
                            neighbors[code] = [node]
                        mem.append(node)
    temp_mem = list(dict.fromkeys(mem))
    mem = temp_mem + path
    return neighbors, mem


def node_json_making(path, output="sp_nodes.txt"):
    print("making nodes file...")
    with open(output, 'w') as filout:
        for i in path:
            i = str(i)
            filout.write(f'{{ data: {{ id: \'{i}\', size: 200, name: \'{i}\' }}, classes: [] }},\n')


def edge_json_making(path, edges, neighbors=None, output="sp_edges.txt"):
    print("making edges file...")
    to_be_wrote = []
    print(" handling shortest path...")
    for i in range(1, len(path)):
        for j in edges:
            listed_edge = j.split("'")
            if (listed_edge[1] == path[i] and listed_edge[3] == path[i - 1]) or (
                    listed_edge[3] == path[i] and listed_edge[1] == path[i - 1]):
                to_be_wrote.append(j)
    if neighbors is not None:
        print(" handling neighbors...")
        for k in edges:
            listed_edge = k.split("'")
            if listed_edge[1] in neighbors.keys():
                if listed_edge[3] in neighbors[listed_edge[1]]:
                    to_be_wrote.append(k)
            if listed_edge[3] in neighbors.keys():
                if listed_edge[1] in neighbors[listed_edge[3]]:
                    to_be_wrote.append(k)
    json = list(dict.fromkeys(to_be_wrote))
    with open(output, 'w') as filout:
        for line in json:
            filout.write(line)


def shortpath_main():
    print("program launched")
    nodes = parse_nodes(args.network[0])
    unicodes = uni_listing(args.codes)
    check = check_nodes(nodes, unicodes)
    if check == 1:
        print("At least one of the proteins is not present")
        exit()
    edges, edges_intel = parse_edges(args.network[1])
    shorter_path = find_path(edges, args.codes[0], args.codes[1], nodes)
    # print(shorter_path)
    if args.neighbors:
        if len(shorter_path) == 1:
            neighbors, true_nodes = neighbor_finder(shorter_path, edges)
            if args.output:
                edge_json_making(shorter_path, edges_intel, neighbors, args.output[1])
                node_json_making(true_nodes, args.output[0])
            else:
                edge_json_making(shorter_path, edges_intel, neighbors)
                node_json_making(true_nodes)
        else:
            for path in shorter_path:
                neighbors, true_nodes = neighbor_finder(path, edges)
                if args.output:
                    edge_json_making(path, edges_intel, neighbors, f"{args.output[1]}_{shorter_path.index(path)}.txt")
                    node_json_making(true_nodes, f"{args.output[0]}_{shorter_path.index(path)}")
                else:
                    edge_json_making(path, edges_intel, neighbors, f"sp_edges_{shorter_path.index(path)}.txt")
                    node_json_making(true_nodes, f"sp_nodes_{shorter_path.index(path)}.txt")
    else:
        if len(shorter_path) == 1:
            if args.output:
                edge_json_making(shorter_path, edges_intel, None, args.output[1])
                node_json_making(shorter_path, args.output[0])
            else:
                edge_json_making(shorter_path, edges_intel)
                node_json_making(shorter_path)
        else:
            for path in shorter_path:
                if args.output:
                    edge_json_making(path, edges_intel, None, f"{args.output[1]}_{shorter_path.index(path)}.txt")
                    node_json_making(path, f"{args.output[0]}_{shorter_path.index(path)}")
                else:
                    edge_json_making(path, edges_intel, None, f"sp_edges_{shorter_path.index(path)}.txt")
                    node_json_making(path, f"sp_nodes_{shorter_path.index(path)}.txt")

    print("program done")


shortpath_main()
