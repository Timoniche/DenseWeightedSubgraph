import unittest
import itertools

from Goldberg import Find_Densest_Subgraph, Find_Density
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density
from utils.graph_generator import gen_random_graph


class TestCase(unittest.TestCase):
    def test_not_weighted(self):
        number_of_nodes = 15
        filepath = "random.txt"
        number_of_edges = gen_random_graph(n=8, filepath=filepath, weighted=False)

        answer_pred = Find_Densest_Subgraph(number_of_nodes, number_of_edges, filepath)
        dens_pred = Find_Density(answer_pred, filepath)

        all_subgraphs = list(itertools.product([0, 1], repeat=number_of_nodes))
        dens_ans = -1
        for subgr in all_subgraphs:
            cur_subgr_indices = []
            for i in range(number_of_nodes):
                if subgr[i] == 0:
                    cur_subgr_indices.append(i)
            if cur_subgr_indices != []:
                dens_ans = max(dens_ans, Find_Density(cur_subgr_indices, filepath))

        self.assertEqual(dens_pred, dens_ans)

    def test_weighted(self):
        number_of_nodes = 15
        filepath = "random_weighted.txt"
        number_of_edges = gen_random_graph(n=8, filepath=filepath, weighted=True)

        answer_pred = WFind_Densest_Subgraph(number_of_nodes, number_of_edges, filepath)
        dens_pred = WFind_Density(answer_pred, filepath)

        all_subgraphs = list(itertools.product([0, 1], repeat=number_of_nodes))
        dens_ans = -1
        for subgr in all_subgraphs:
            cur_subgr_indices = []
            for i in range(number_of_nodes):
                if subgr[i] == 0:
                    cur_subgr_indices.append(i)
            if cur_subgr_indices != []:
                dens_ans = max(dens_ans, WFind_Density(cur_subgr_indices, filepath))

        self.assertEqual(dens_pred, dens_ans)


if __name__ == '__main__':
    unittest.main()
