import random


# returns edges cnt
def gen_random_graph(n, filepath, weighted):
    file = open(filepath, "w")
    adj = []
    for i in range(n):
        i_set = set()
        runs_cnt = random.randint(0, n - 1)
        for j in range(runs_cnt):
            v_to = random.randint(0, n - 1)
            if v_to != i:
                i_set.add(v_to)
        adj.append(i_set)
    edges = set()
    for i in range(n):
        for v_to in adj[i]:
            edges.add((min(i, v_to), max(i, v_to)))
    for e in edges:
        if not weighted:
            file.write(f'{e[0]} {e[1]}\n')
        else:
            rand_dist = random.uniform(1, 5)
            file.write(f'{e[0]} {e[1]} {rand_dist}\n')

    file.close()
    return len(edges)
