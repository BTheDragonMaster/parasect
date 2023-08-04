from argparse import ArgumentParser
from dataclasses import dataclass
import typing as ty
import os

import torch
from torch import Tensor
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.optim import Optimizer, lr_scheduler, Adam
from tqdm import tqdm
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
import matplotlib.pyplot as plt

from paras.scripts.data_processing.structure_processing.make_voxel_dataset import VoxelDataset
import paras.data.structure_data

SUBSTRATE_PATH = os.path.join(os.path.dirname(paras.data.structure_data.__file__), "included_substrates.txt")


class Cnn3D(torch.nn.Module):
    pass


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-train', type=str, required=True, help="Domain list of training data points")
    parser.add_argument('-test', type=str, required=True, help="Domain list of test data points")
    parser.add_argument('-i', type=str, required=True, help="Directory to saved voxel features directory")
    parser.add_argument('-s', type=str, required=True, help="Directory to specificities file")
    parser.add_argument('-e', type=int, default=35, help="Number of epochs")
    parser.add_argument('-b', type=int, default=10, help="Batch size")
    parser.add_argument('-o', type=str, required=True, help="Output dir")
    parser.add_argument('-m', type=str, default=None, help="Directory to pre-trained model")
    parser.add_argument('-l', type=float, default=0.001, help="Start learning rate")
    parser.add_argument('-n', type=str, default=SUBSTRATE_PATH, help="Path to file containing included substrates")

    args = parser.parse_args()
    return args


@dataclass
class SubstrateMetrics:
    substrate: str

    tp_count: int = 0
    fp_count: int = 0
    fn_count: int = 0
    tn_count: int = 0

    def __post_init__(self):
        self.mispredictions = []

    def process_prediction(self, pred_y, true_y):
        if pred_y > 0.5 and true_y > 0.5:
            self.tp_count += 1
        elif pred_y <= 0.5 < true_y:
            self.fn_count += 1
        elif true_y <= 0.5 < pred_y:
            self.fp_count += 1
        else:
            self.tn_count += 1

    def count_mispredictions(self):
        counts = {}
        for misprediction in self.mispredictions:
            if misprediction not in counts:
                counts[misprediction] = 0
            counts[misprediction] += 1

        return counts

    def print_mispredictions(self):
        counts = self.count_mispredictions()
        print("Mispredictions")
        sorted_counts = []
        for misprediction, count in counts.items():
            sorted_counts.append((misprediction, count))
        sorted_counts.sort(key=lambda x: x[1], reverse=True)

        for misprediction, count in sorted_counts:
            print(f"{misprediction}: {count}")

    def print(self):
        print(f"Substrate: {self.substrate}")
        print(f"TP: {self.tp_count}")
        print(f"FP: {self.fp_count}")
        print(f"TN: {self.tn_count}")
        print(f"FN: {self.fn_count}\n")
        print('\n')
        self.print_mispredictions()
        print('\n')


def write_aa_metrics(out_folder, epoch, metrics, specificities):
    out_file = os.path.join(out_folder, f"epoch_{epoch}.txt")
    out_file_2 = os.path.join(out_folder, f"epoch_{epoch}_cm.txt")
    with open(out_file, 'w') as out:
        out.write("amino_acid\tTP\tFP\tTN\tFN\n")
        for aa in specificities:
            aa_metrics = metrics[aa]
            out.write(f"{aa}\t{aa_metrics.tp_count}\t{aa_metrics.fp_count}\t{aa_metrics.tn_count}\t{aa_metrics.fn_count}\n")
    with open(out_file_2, 'w') as out:
        out.write('Predicted\\Actual\t')
        out.write('\t'.join(specificities))
        out.write("\n")
        for aa in specificities:
            out.write(f"{aa}")
            aa_metrics = metrics[aa]
            misprediction_counts = aa_metrics.count_mispredictions()
            for aa_2 in specificities:
                if aa_2 not in misprediction_counts:
                    if aa_2 == aa:
                        out.write(f"\t{aa_metrics.tp_count}")
                    else:
                        out.write(f"\t{0}")
                else:
                    out.write(f"\t{misprediction_counts[aa_2]}")
            out.write("\n")


def plot_aa_metrics(out_folder, metrics, specificities):
    for specificity in specificities:
        out_file = os.path.join(out_folder, f"{specificity}.svg")
        x = []
        tp = []
        fp = []
        tn = []
        fn = []

        for epoch, aa_to_metrics in metrics.items():
            x.append(epoch)

            tp.append(aa_to_metrics[specificity].tp_count)
            fp.append(aa_to_metrics[specificity].fp_count)
            tn.append(aa_to_metrics[specificity].tn_count)
            fn.append(aa_to_metrics[specificity].fn_count)

        plt.plot(x, tp, label="TP")
        plt.plot(x, fp, label="FP")
        plt.plot(x, tn, label="TN")
        plt.plot(x, fn, label="FN")
        plt.legend()
        plt.savefig(out_file)
        plt.clf()


def write_metrics(out_file, epochs, train_losses, test_losses, accuracies, learning_rates, precisions, recalls, class_precisions, class_recalls, train_accuracies):
    with open(out_file, 'w') as out:
        out.write("epoch\ttrain_loss\ttest_loss\taccuracy\tlearning_rate\tprecision\trecall\tclass_precision\tclass_recall\ttrain_accuracy\n")
        for i, epoch in enumerate(epochs):
            out.write(f"{epoch}\t{train_losses[i]}\t{test_losses[i]}\t{accuracies[i]}\t{learning_rates[i]}\t{precisions[i]}\t{recalls[i]}\t{class_precisions[i]}\t{class_recalls[i]}\t{train_accuracies[i]}\n")


def plot_metrics(out_plots, epochs, train_losses, test_losses, accuracies, precisions, recalls, class_precisions, class_recalls, train_accuracies):
    precision_out = os.path.join(out_plots, "precision.png")
    recall_out = os.path.join(out_plots, "recall.png")
    class_precision_out = os.path.join(out_plots, "class_precision.png")
    class_recall_out = os.path.join(out_plots, "class_recall.png")
    accuracy_out = os.path.join(out_plots, "accuracy.png")
    loss_out = os.path.join(out_plots, "loss.png")

    plt.plot(epochs, precisions)
    plt.savefig(precision_out)
    plt.clf()

    plt.plot(epochs, recalls)
    plt.savefig(recall_out)
    plt.clf()

    plt.plot(epochs, class_precisions)
    plt.savefig(class_precision_out)
    plt.clf()

    plt.plot(epochs, class_recalls)
    plt.savefig(class_recall_out)
    plt.clf()

    plt.plot(epochs, train_accuracies, label="Train")
    plt.plot(epochs, accuracies, label="Test")
    plt.legend()
    plt.savefig(accuracy_out)
    plt.clf()

    plt.plot(epochs, train_losses, label="Train")
    plt.plot(epochs, test_losses, label="Test")
    plt.legend()
    plt.savefig(loss_out)
    plt.clf()


class Convolution3DSmall(torch.nn.Module):
    def __init__(
            self,
            n_classes: int,
            in_channels: int = 7,
            out_channels: int = 16,
            kernel_size: int = 3,
            pool_size: int = 2

    ) -> None:
        super().__init__()
        self.pool_size = pool_size
        self.kernel_size = tuple([kernel_size] * 3)
        self.conv1 = self.make_conv_layer_pool(in_channels, out_channels)
        self.conv2 = self.make_conv_layer_pool(out_channels, out_channels * 2)
        self.fc1 = self.make_fc_layer(2**3 * out_channels * 2, 128, dropout=True)
        self.fc2 = self.make_fc_layer(128, n_classes)

    def make_conv_layer(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU()
        )

        return conv_layer

    def make_conv_layer_pool(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU(),
            nn.AvgPool3d(self.pool_size, stride=self.pool_size, padding=0, count_include_pad=False)
        )

        return conv_layer

    def make_fc_layer(self, in_channels, out_channels, dropout=False):
        if dropout:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels),
                nn.Dropout(p=0.15)

            )
        else:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels)
            )
        return fc_layer

    def forward(self, x: Tensor) -> torch.Tensor:
        out = self.conv1(x)
        out = self.conv2(out)
        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        out = self.fc2(out)

        return out


class Convolution3D(torch.nn.Module):
    def __init__(
            self,
            n_classes: int,
            input_size: int,
            in_channels: int = 7,
            conv_layers_out: ty.Tuple = (32, 32, 32),
            conv_pool: ty.Tuple = (2, 0, 0),
            conv_padding: ty.Tuple = (1, 1, 0),
            conv_stride: ty.Tuple = (1, 1, 1),
            hidden_fc_layers_out: ty.Tuple = (512, 256, 128),
            fc_dropout: ty.Tuple = (0.15, 0.15, 0.15),
            kernel_size: ty.Tuple = (3, 3, 3),

    ) -> None:
        super().__init__()
        assert len(conv_layers_out) == len(conv_pool) == len(conv_padding) == len(conv_stride)
        assert len(hidden_fc_layers_out) == len(fc_dropout)

        self.conv_layers = []

        nr_channels = in_channels
        box_size = input_size
        out_channels = in_channels

        for i, conv_layer_out in enumerate(conv_layers_out):
            conv_layer = self.make_conv_layer(nr_channels, conv_layer_out,
                                              padding=conv_padding[i],
                                              stride=conv_stride[i],
                                              pool=conv_pool[i], kernel_size=kernel_size[i])
            box_size = self.calculate_box_size(box_size, conv_stride[i], conv_pool[i], conv_padding[i], kernel_size[i])
            self.conv_layers.append(conv_layer)
            out_channels = conv_layer_out

        out_channels = box_size ** 3 * out_channels

        self.hidden_layers = []

        for i, hidden_layer_out in enumerate(hidden_fc_layers_out):
            dropout = fc_dropout[i]
            hidden_layer = self.make_fc_layer(out_channels, hidden_layer_out, dropout=dropout)
            out_channels = hidden_layer_out
            self.hidden_layers.append(hidden_layer)

        self.output_layer = self.make_fc_layer(out_channels, n_classes, dropout=0.0)

    def calculate_box_size(self, previous_box, stride, pool, padding, kernel_size):
        box_size = ((previous_box - kernel_size + 2 * padding) / stride) + 1
        if pool:
            box_size = ((box_size - self.pool_size) / self.pool_size) + 1

        return box_size

    def make_conv_layer(self, in_channels, out_channels, padding=0, stride=1, pool=0, kernel_size=3):
        conv_layer = nn.Conv3d(in_channels, out_channels, kernel_size=tuple([kernel_size] * 3), padding=padding,
                               stride=tuple([stride] * 3))
        activation_layer = nn.ReLU()

        if not pool:
            full_conv_layer = nn.Sequential(conv_layer, activation_layer)

        else:
            pool_layer = nn.AvgPool3d(pool, stride=pool, padding=0, count_include_pad=False)
            full_conv_layer = nn.Sequential(conv_layer, activation_layer, pool_layer)

        return full_conv_layer

    @staticmethod
    def make_fc_layer(in_channels, out_channels, dropout=0.15):
        fc_layer = nn.Sequential(nn.Linear(in_channels, out_channels),
                                 nn.ReLU(),
                                 nn.BatchNorm1d(out_channels),
                                 nn.Dropout(p=dropout))
        return fc_layer

    def forward(self, x: Tensor) -> torch.Tensor:
        out = x
        for conv_layer in self.conv_layers:
            out = conv_layer(out)

        out = out.view(out.size(0), -1)
        for hidden_layer in self.hidden_layers:
            out = hidden_layer(out)

        out = self.output_layer(out)

        return out


class Convolution3DOld(torch.nn.Module):
    def __init__(
            self,
            n_classes: int,
            in_channels: int = 7,
            out_channels: int = 32,
            kernel_size: int = 3,
            pool_size: int = 2

    ) -> None:
        super().__init__()
        self.pool_size = pool_size
        self.kernel_size = tuple([kernel_size] * 3)
        self.conv1 = self.make_conv_layer(in_channels, out_channels, padding=1)
        self.conv2 = self.make_conv_layer(out_channels, out_channels, padding=1)
        self.conv3 = self.make_conv_layer(out_channels, out_channels, padding=1)
        self.fc1 = self.make_fc_layer(17**3 * out_channels, 512)
        self.fc2 = self.make_fc_layer(512, 256)
        self.fc3 = self.make_fc_layer(256, 128)
        self.fc4 = self.make_fc_layer(128, n_classes, dropout=False)

    def make_conv_layer(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU()
        )

        return conv_layer

    def make_conv_layer_pool(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU(),
            nn.AvgPool3d(self.pool_size, stride=self.pool_size, padding=0, count_include_pad=False)
        )

        return conv_layer

    def make_fc_layer(self, in_channels, out_channels, dropout=True):
        if dropout:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels),
                nn.Dropout(p=0.15)

            )
        else:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels))

        return fc_layer

    def forward(self, x: Tensor) -> torch.Tensor:
        out = self.conv1(x)
        out = self.conv2(out)
        out = self.conv3(out)

        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        out = self.fc2(out)
        out = self.fc3(out)
        out = self.fc4(out)

        return out


class Convolution3DSeq(torch.nn.Module):
    def __init__(
            self,
            n_classes: int,
            in_channels: int = 7,
            out_channels: int = 32,
            kernel_size: int = 3,
            pool_size: int = 2

    ) -> None:
        super().__init__()
        self.pool_size = pool_size
        self.kernel_size = tuple([kernel_size] * 3)
        self.conv1 = self.make_conv_layer(in_channels, out_channels, padding=1)
        self.conv2 = self.make_conv_layer(out_channels, out_channels, padding=1)
        self.conv3 = self.make_conv_layer(out_channels, out_channels, padding=1)
        self.fc1 = self.make_fc_layer(17**3 * out_channels, 512)
        self.fc2 = self.make_fc_layer(512, 256)
        self.fc3 = self.make_fc_layer(256, 128)
        self.fc4 = self.make_fc_layer(128, n_classes, dropout=False)

    def make_conv_layer(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU()
        )

        return conv_layer

    def make_conv_layer_pool(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU(),
            nn.AvgPool3d(self.pool_size, stride=self.pool_size, padding=0, count_include_pad=False)
        )

        return conv_layer

    def make_fc_layer(self, in_channels, out_channels, dropout=True):
        if dropout:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels),
                nn.Dropout(p=0.15)

            )
        else:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels))

        return fc_layer

    def forward(self, x: Tensor) -> torch.Tensor:
        out = self.conv1(x)
        out = self.conv2(out)
        out = self.conv3(out)

        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        out = self.fc2(out)
        out = self.fc3(out)
        out = self.fc4(out)

        return out


class Convolution3DLarge(torch.nn.Module):
    def __init__(
            self,
            n_classes: int,
            in_channels: int = 7,
            out_channels: int = 32,
            kernel_size: int = 3,
            pool_size: int = 2

    ) -> None:
        super().__init__()
        self.pool_size = pool_size
        self.kernel_size = tuple([kernel_size] * 3)
        self.conv1 = self.make_conv_layer_pool(in_channels, out_channels, padding=1)
        self.conv2 = self.make_conv_layer(out_channels, out_channels * 2, padding=1)
        self.conv3 = self.make_conv_layer(out_channels * 2, out_channels * 2, padding=1)
        self.fc1 = self.make_fc_layer(8**3 * out_channels * 2, 16384)
        self.fc2 = self.make_fc_layer(16384, 4096)
        self.fc3 = self.make_fc_layer(4096, 1024)
        self.fc4 = self.make_fc_layer(1024, 256, dropout=False)
        self.fc5 = self.make_fc_layer(256, 128, dropout=False)
        self.fc6 = self.make_fc_layer(128, n_classes, dropout=False)

    def make_conv_layer(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU()
        )

        return conv_layer

    def make_conv_layer_pool(self, in_channels, out_channels, padding=0):
        conv_layer = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=self.kernel_size, padding=padding),
            nn.ReLU(),
            nn.AvgPool3d(self.pool_size, stride=self.pool_size, padding=0, count_include_pad=False)
        )

        return conv_layer

    def make_fc_layer(self, in_channels, out_channels, dropout=True):
        if dropout:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels),
                nn.Dropout(p=0.15)

            )
        else:
            fc_layer = nn.Sequential(
                nn.Linear(in_channels, out_channels),
                nn.ReLU(),
                nn.BatchNorm1d(out_channels))

        return fc_layer

    def forward(self, x: Tensor) -> torch.Tensor:
        out = self.conv1(x)
        out = self.conv2(out)
        out = self.conv3(out)
        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        out = self.fc2(out)
        out = self.fc3(out)
        out = self.fc4(out)
        out = self.fc5(out)
        out = self.fc6(out)

        return out


def train_loop(
        model: nn.Module,
        loader: DataLoader,
        optimizer: Optimizer,
        criterion: ty.Callable
) -> ty.Tuple[np.ndarray, np.ndarray]:
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
    loader = tqdm(loader, leave=False)
    batch_losses = []
    accuracies = []
    for protein, specificities in loader:
        optimizer.zero_grad()
        out = model(protein)
        result = calculate_metrics(out.tolist(), specificities.tolist())
        loss_batch = criterion(out, specificities.type(torch.float))
        batch_loss_value = loss_batch.item()
        loss_batch.backward()
        optimizer.step()
        batch_losses.append(batch_loss_value)
        accuracies.append(result['accuracy'])
    loss_value = np.mean(batch_losses)
    mean_accuracy = np.mean(np.array(accuracies))

    return loss_value, mean_accuracy


def eval_loop(
        model: nn.Module,
        loader: DataLoader,
        criterion: ty.Callable,
        epoch: int,
        metrics: ty.Dict[int, ty.Dict[str, SubstrateMetrics]],
        specificity_list: ty.List[str]
) -> ty.Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluation loop for Shapenet.

    Arguments
    ---------
    model (nn.Module): Neural net model.
    loader (DataLoader): Data loader.
    criterion (ty.Callable): Loss function.

    Returns
    -------
    ty.Tuple[float, float]: Test loss and accuracy.
    """
    model.eval()
    loader = tqdm(loader, leave=False)
    batch_losses = []
    batch_nr = 0
    accuracies = []
    recalls = []
    class_recalls = []
    precisions = []
    class_precisions = []

    for protein, specificities in loader:
        batch_nr += 1
        with torch.no_grad():
            prediction = model(protein)
            loss_batch = criterion(prediction, specificities.type(torch.float))
            batch_loss_value = loss_batch.item()
            batch_losses.append(batch_loss_value)
            result = calculate_metrics(prediction.tolist(), specificities.tolist())
            store_aa_metrics(metrics, prediction, specificities, specificity_list, epoch)

            recalls.append(result['micro/recall'])
            precisions.append(result['micro/precision'])
            class_recalls.append(result['macro/recall'])
            class_precisions.append(result['macro/precision'])
            accuracies.append(result['accuracy'])

    loss_value = np.mean(batch_losses)
    mean_accuracy = np.mean(np.array(accuracies))
    mean_recall = np.mean(np.array(recalls))
    mean_precision = np.mean(np.array(precisions))
    mean_class_recall = np.mean(np.array(class_recalls))
    mean_class_precision = np.mean(np.array(class_precisions))

    return loss_value, mean_accuracy, mean_recall, mean_precision, mean_class_recall, mean_class_precision


def save_model(
        epoch: int,
        model: nn.Module,
        optimizer: Optimizer,
        save_dir: str,
) -> None:
    """
    Save model.

    Arguments
    ---------
    epoch (int): Epoch.
    model (nn.Module): Neural net model.
    optimizer (Optimizer): Optimizer.
    save_dir (str): Save directory.
    """
    state = {
        "epoch": epoch,
        "state_dict": model.state_dict(),
        "optimizer": optimizer.state_dict()
    }
    torch.save(state, os.path.join(save_dir, f"model_epoch_{epoch + 1}.pth.tar"))


def load_checkpoint(model, optimizer, filename):
    # Note: Input model & optimizer should be pre-defined.  This routine only updates their states.
    start_epoch = 0
    if os.path.isfile(filename):
        print("=> loading checkpoint '{}'".format(filename))
        checkpoint = torch.load(filename)
        start_epoch = checkpoint['epoch']
        model.load_state_dict(checkpoint['state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        print(f"=> loaded checkpoint '{filename}' (epoch {start_epoch})")
    else:
        print(f"=> no checkpoint found at '{filename}'")

    return model, optimizer, start_epoch


def main() -> None:
    """
    Driver code.
    """
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    train = VoxelDataset(args.i, args.train, args.s, SUBSTRATE_PATH)
    test = VoxelDataset(args.i, args.test, args.s, SUBSTRATE_PATH)

    if len(train) % args.b == 1:

        train_loader = DataLoader(train, shuffle=True, batch_size=args.b, drop_last=True)
    else:
        train_loader = DataLoader(train, shuffle=True, batch_size=args.b)

    if len(test) % args.b == 1:
        test_loader = DataLoader(test, shuffle=True, batch_size=args.b, drop_last=True)
    else:
        test_loader = DataLoader(test, shuffle=True, batch_size=args.b)

    model = Convolution3D(len(train.all_specificities))

    weights = torch.FloatTensor(train.get_weights())

    lr = args.l
    criterion = nn.BCEWithLogitsLoss(pos_weight=weights)
    optimizer = Adam(params=model.parameters(), lr=lr)

    start_epoch = 0

    if args.m:
        model, optimizer, start_epoch = load_checkpoint(model, optimizer, args.m)

    scheduler = lr_scheduler.ReduceLROnPlateau(
        optimizer=optimizer,
        mode='min',
        factor=0.5,
        patience=5,
        min_lr=0.00001
    )

    epochs, train_losses, test_losses, accuracies, learning_rates, precisions, recalls, class_precisions, class_recalls, train_accuracies = \
        [], [], [], [], [], [], [], [], [], []

    num_epochs = args.e
    aa_metrics_folder = os.path.join(args.o, "amino_acid_metrics")
    if not os.path.exists(aa_metrics_folder):
        os.mkdir(aa_metrics_folder)
    aa_metrics = {}
    for epoch in range(num_epochs):
        aa_metrics[epoch + 1 + start_epoch] = {}
        for spec in train.all_specificities:
            aa_metrics[epoch + 1 + start_epoch][spec] = SubstrateMetrics(spec)

    for epoch in range(num_epochs):

        lr = scheduler.optimizer.param_groups[0]['lr']
        train_loss, train_accuracy = train_loop(
            model=model,
            loader=train_loader,
            optimizer=optimizer,
            criterion=criterion
        )
        test_loss, accuracy, recall, precision, class_recall, class_precision = eval_loop(
            model=model,
            loader=test_loader,
            criterion=criterion,
            epoch=epoch + start_epoch + 1,
            metrics=aa_metrics,
            specificity_list=train.all_specificities
        )
        scheduler.step(test_loss)
        for aa, metrics in aa_metrics[epoch + start_epoch + 1].items():
            metrics.print()

        epochs.append(epoch + 1 + start_epoch)
        train_losses.append(train_loss)
        test_losses.append(test_loss)
        learning_rates.append(lr)
        accuracies.append(accuracy)
        recalls.append(recall)
        class_recalls.append(class_recall)
        precisions.append(precision)
        class_precisions.append(class_precision)
        train_accuracies.append(train_accuracy)

        write_aa_metrics(aa_metrics_folder, epoch + start_epoch + 1, aa_metrics[epoch + start_epoch + 1], train.all_specificities)

        print(
            f"Epoch: {epoch + 1 + start_epoch}/{num_epochs + start_epoch}; "
            f"LR: {lr}; "
            f"Train loss: ", train_loss,
            f"Test loss: ", test_loss,
            f"Accuracy: ", accuracy,
            f"Recall: ", recall,
            f"Precision: ", precision,
        )

        if (epoch + 1) % 10 == 0:
            save_model(epoch + start_epoch, model, optimizer, args.o)

    out_plots = os.path.join(args.o, "plots")

    if not os.path.exists(out_plots):
        os.mkdir(out_plots)
    metric_plots = os.path.join(out_plots, "amino_acid_metrics")
    if not os.path.exists(metric_plots):
        os.mkdir(metric_plots)

    plot_aa_metrics(metric_plots, aa_metrics, train.all_specificities)

    metrics_out = os.path.join(args.o, "metrics.txt")
    write_metrics(metrics_out, epochs, train_losses, test_losses, accuracies, learning_rates, precisions, recalls,
                  class_precisions, class_recalls, train_accuracies)

    plot_metrics(out_plots, epochs, train_losses, test_losses, accuracies, precisions, recalls,
                 class_precisions, class_recalls, train_accuracies)

    save_model(num_epochs, model, optimizer, args.o)


def calculate_metrics(pred, target, threshold=0.5):
    pred = np.array(pred)
    target = np.array(target)
    pred = np.array(pred > threshold, dtype=float)

    return {'micro/precision': precision_score(y_true=target, y_pred=pred, average='micro', zero_division=0),
            'micro/recall': recall_score(y_true=target, y_pred=pred, average='micro', zero_division=0),
            'micro/f1': f1_score(y_true=target, y_pred=pred, average='micro', zero_division=0),
            'macro/precision': precision_score(y_true=target, y_pred=pred, average='macro', zero_division=0),
            'macro/recall': recall_score(y_true=target, y_pred=pred, average='macro', zero_division=0),
            'macro/f1': f1_score(y_true=target, y_pred=pred, average='macro', zero_division=0),
            'samples/precision': precision_score(y_true=target, y_pred=pred, average='samples', zero_division=0),
            'samples/recall': recall_score(y_true=target, y_pred=pred, average='samples', zero_division=0),
            'samples/f1': f1_score(y_true=target, y_pred=pred, average='samples', zero_division=0),
            'accuracy': accuracy_score(y_true=target, y_pred=pred)
            }


def store_aa_metrics(metrics, prediction, specificities, specificity_list, epoch):

    for i, pred in enumerate(prediction):
        spec = specificities[i]
        correct_indices = []

        for j, true_y in enumerate(spec):
            pred_y = pred[j]
            specificity = specificity_list[j]
            metrics[epoch][specificity].process_prediction(pred_y, true_y)

            if true_y > 0.5:
                correct_indices.append(j)

        for k, pred_y in enumerate(pred):
            if pred_y > 0.5 and k not in correct_indices:
                for idx in correct_indices:
                    metrics[epoch][specificity_list[k]].mispredictions.append(specificity_list[idx])


if __name__ == "__main__":
    main()
