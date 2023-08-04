import torch
import torch.nn.functional as F
from torch.nn import Sigmoid
from torch_geometric.nn import GCNConv, Linear
from torch_geometric.nn.aggr import Set2Set
from torch_geometric.data import Data, Batch
from torch_geometric.loader import DataLoader



class NeuralNet(torch.nn.Module):
    def __init__(self):
        super().__init__()

        self.activation = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(0.3)

    def forward(self, enzyme, molecule):
        pass


class PairData(Data):
    def __init__(self, edge_index_s=None, x_s=None, edge_index_t=None, x_t=None):
        super().__init__()
        self.edge_index_s = edge_index_s
        self.x_s = x_s
        self.edge_index_t = edge_index_t
        self.x_t = x_t

    def __inc__(self, key, value, *args, **kwargs):
        if key == 'edge_index_s':
            return self.x_s.size(0)
        if key == 'edge_index_t':
            return self.x_t.size(0)
        else:
            return super().__inc__(key, value, *args, **kwargs)


class HeteroPairData(Data):
    def __init__(self, edge_index_s=None, x_s=None, x_t=None):
        super().__init__()
        self.edge_index_s = edge_index_s
        self.x_s = x_s
        self.x_t = x_t

    def __inc__(self, key, value, *args, **kwargs):
        if key == 'edge_index_s':
            return self.x_s.size(0)
        else:
            return super().__inc__(key, value, *args, **kwargs)


if __name__ == "__main__":
    edge_index_s = torch.tensor([
        [0, 0, 0, 0],
        [1, 2, 3, 4],
    ])
    x_s = torch.randn(5, 17)  # 5 nodes.
    edge_index_t = torch.tensor([
        [0, 0, 0],
        [1, 2, 3],
    ])
    x_t = torch.randn(7, 17, 17, 17)  # 4 nodes.

    data = HeteroPairData(edge_index_s, x_s, x_t)
    data_list = [data, data, data, data, data, data]
    # loader = DataLoader(data_list, batch_size=2)
    # batch = next(iter(loader))
    #
    # print(batch)
    #
    # print(batch.edge_index_s)
    #
    # print(batch.edge_index_t)

    loader = DataLoader(data_list, batch_size=3, follow_batch=['x_s', 'x_t'])
    for i, batch in enumerate(loader):
        print(i)

        # batch = next(iter(loader))
        print(batch)
        print(batch.x_s_batch)
        print(batch.x_t_batch)
        batch = next(iter(loader))
        print(batch)
        print(batch.x_s_batch)
        print(batch.x_t_batch)
        print("HERE")
        print(batch.edge_index_s)

        print(batch.x_s.size())
        print(batch.x_t.size())
        # set2set = Set2Set(5, 1)
        # print(set2set(batch.x_t, batch.x_t_batch))
        # lin = Linear(10, 1)
        # output = lin(set2set(batch.x_t, batch.x_t_batch))
        # sigmoid = Sigmoid()
        # print(sigmoid(output))

