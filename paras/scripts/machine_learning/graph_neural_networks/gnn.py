import torch
import torch.nn.functional as F
import torch.nn as nn
from torch_geometric.nn import NNConv, GCNConv
from torch_geometric.nn.aggr import Set2Set
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from tqdm import tqdm


class MPNN(torch.nn.Module):
    def __init__(
            self,
            input_size_nodes: int,
            input_size_edges: int,
            hidden_size: int,
            embed_size: int,
            message_passing_steps: int = 3
    ) -> None:
        super().__init__()
        self.message_passing_steps = message_passing_steps
        self.project = nn.Linear(input_size_nodes, hidden_size)
        edge_net = nn.Sequential(
            nn.Linear(input_size_edges, 2 * hidden_size),
            nn.ReLU(),
            nn.Linear(2 * hidden_size, hidden_size * hidden_size)
        )
        self.conv = NNConv(hidden_size, hidden_size, edge_net, aggr="mean")  # Maps edge features onto nodes
        self.gru = nn.GRU(hidden_size, hidden_size)
        self.set2set = Set2Set(hidden_size, processing_steps=3)
        self.embed = nn.Linear(2 * hidden_size, embed_size)

    def forward(self, data: Data) -> torch.Tensor:
        out = F.relu(self.project(data.x))
        h = out.unsqueeze(0)  # adds a dimension; output without edges; needed for GRU
        for _ in range(self.message_passing_steps):
            m = F.relu(self.conv(out, data.edge_index, data.edge_attr))  # Output with edges; collects important edges (for output) (salt)
            out, h = self.gru(m.unsqueeze(0), h)  # Gated Recurrent Unit -> Max Welling paper
            out = out.squeeze(0)
        out = self.set2set(out, data.batch)
        return F.relu(self.embed(out))


class Model(torch.nn.Module):
    def __init__(
            self,

            # Parameters for substrate embedder:
            input_size_nodes_substrate: int,
            input_size_edges_substrate: int,
            hidden_size_substrate: int,
            embed_size_substrate: int,
            message_passing_steps_substrate: int,

            # Parameters for enzyme pocket embedder.
            input_size_nodes_enzyme: int,  # = 14
            input_size_edges_enzyme: int,  # = 1
            hidden_size_enzyme: int,
            embed_size_enzyme: int,
            message_passing_steps_enzyme: int,
    ) -> None:
        super().__init__()
        self.substrate_embedder = MPNN(
            input_size_nodes=input_size_nodes_substrate,
            input_size_edges=input_size_edges_substrate,
            hidden_size=hidden_size_substrate,
            embed_size=embed_size_substrate,
            message_passing_steps=message_passing_steps_substrate
        )
        self.enzyme_embedder = MPNN(
            input_size_nodes=input_size_nodes_enzyme,
            input_size_edges=input_size_edges_enzyme,
            hidden_size=hidden_size_enzyme,
            embed_size=embed_size_enzyme,
            message_passing_steps=message_passing_steps_enzyme
        )
        total_embed_size = embed_size_substrate + embed_size_enzyme
        self.predict = nn.Sequential(
            nn.Linear(total_embed_size, int(total_embed_size / 2)),
            nn.ReLU(),
            nn.Linear(int(total_embed_size / 2), 1)
        )

    def forward(self, substrate: Data, enzyme: Data) -> torch.Tensor:
        out_substrate = self.substrate_embedder(substrate)
        out_enzyme = self.enzyme_embedder(enzyme)
        out = torch.cat([out_substrate, out_enzyme], dim=1)
        return self.predict(out)


def train_loop(
    model: nn.Module,
    loader: DataLoader,
    optimizer: Optimizer,
    criterion: ty.Callable
) -> float:
    """
    Train loop for Shapenet.
    Arguments
    ---------
    model (nn.Module): Neural net model.
    loader (DataLoader): Data loader.
    optimizer (Optimizer): Optimizer.
    criterion (ty.Callable): Loss function.
    Returns
    -------
    float: Train loss.
    """
    model.train()
    num_samples = len(loader.dataset)
    loader = tqdm(loader, leave=False)
    acc_loss = 0.0
    for batch in loader:
        optimizer.zero_grad()
        out = model(batch).squeeze()
        loss_batch = criterion(out, batch.y)
        loss_batch.backward()
        optimizer.step()
        acc_loss += loss_batch.item()
    avg_loss = acc_loss / num_samples
    return avg_loss
