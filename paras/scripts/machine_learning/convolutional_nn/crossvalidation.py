from argparse import ArgumentParser
import os

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.optim import lr_scheduler, Adam

from paras.scripts.data_processing.structure_processing.make_voxel_dataset import VoxelDataset
from paras.scripts.machine_learning.convolutional_nn.convolution_3d import train_loop, eval_loop, SUBSTRATE_PATH, \
    SubstrateMetrics, write_metrics, plot_metrics, plot_aa_metrics, save_model, Convolution3D, write_aa_metrics
from paras.scripts.parsers.iterate_over_dir import find_crossval_pairs


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-c', type=str, required=True, help="Directory containing crossvalidation domain lists")
    parser.add_argument('-i', type=str, required=True, help="Directory to saved voxel features directory")
    parser.add_argument('-s', type=str, required=True, help="Directory to specificities file")
    parser.add_argument('-o', type=str, required=True, help="Output dir")

    parser.add_argument('-g', type=int, default=17, help="Input grid size")
    parser.add_argument('-n', type=str, default=SUBSTRATE_PATH, help="Path to file containing included substrates")

    parser.add_argument('-e', type=int, default=100, help="Number of epochs")
    parser.add_argument('-b', type=int, default=10, help="Batch size")
    parser.add_argument('-l', type=float, default=0.001, help="Start learning rate")
    parser.add_argument('-w', action="store_true", help="If toggled, use binary weights")

    parser.add_argument('-nconv', type=int, default=3, help="Number of convolutional layers")
    parser.add_argument('-stride', type=int, nargs='?', default=(1, 1, 1),
                        help="Stride used in each convolutional layer. Ensure the number of values given here are \
                        equal to the nconv parameter.")
    parser.add_argument('-pool', type=int, nargs='?', default=(2, 0, 0),
                        help="Pool size for each convolutional layer, 0 if no pooling layer is used. Ensure the number \
                        of values given here are equal to the nconv parameter.")
    parser.add_argument('-padding', type=int, nargs='?', default=(1, 1, 0),
                        help="Padding used in each convolutional layer. Ensure the number of values given here are \
                        equal to the nconv parameter.")
    parser.add_argument('-convout', type=int, nargs='?', default=(32, 32, 32),
                        help="Output channels for each convolutional layer. Ensure the number of values given here are \
                        equal to the nconv parameter.")
    parser.add_argument('-kernel', type=int, nargs='?', default=(3, 3, 3),
                        help="Kernel sizes for convolutional layer. Ensure the number of values given here are equal \
                        to the nconv parameter.")

    parser.add_argument('-nfc', type=int, default=3,
                        help="Number of hidden fully connected layers.")
    parser.add_argument('-dropout', type=float, nargs='?', default=(0.15, 0.15, 0.15),
                        help="Dropout used in each fully connected layer. Ensure the number of values given here are \
                        equal to the nfc parameter.")
    parser.add_argument('-fcout', type=int, nargs='?', default=(512, 256, 128),
                        help="Output channels for each hidden layer. Ensure the number of values given here are equal \
                        to the nfc parameter.")

    args = parser.parse_args()
    return args





def crossvalidation() -> None:
    """
    Driver code.
    """
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    assert args.nconv == len(args.stride) == len(args.pool) == len(args.padding) == len(args.convout) == \
           len(args.kernel)
    assert args.nfc == len(args.dropout) == len(args.fcout)

    crossval_pairs = find_crossval_pairs(args.c)

    for i, crossval_pair in enumerate(crossval_pairs):
        output_dir = os.path.join(args.o, f"crossval_{i:02d}")

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        train_domain_path, test_domain_path = crossval_pair
        train = VoxelDataset(args.i, train_domain_path, args.s, SUBSTRATE_PATH, grid_size=args.g)
        test = VoxelDataset(args.i, test_domain_path, args.s, SUBSTRATE_PATH, grid_size=args.g)

        if len(train) % args.b == 1:
            train_loader = DataLoader(train, shuffle=True, batch_size=args.b, drop_last=True)
        else:
            train_loader = DataLoader(train, shuffle=True, batch_size=args.b)

        if len(test) % args.b == 1:
            test_loader = DataLoader(test, shuffle=True, batch_size=args.b, drop_last=True)
        else:
            test_loader = DataLoader(test, shuffle=True, batch_size=args.b)

        model = Convolution3D(len(train.all_specificities), args.g,
                              conv_layers_out=args.convout,
                              conv_pool=args.pool,
                              conv_padding=args.padding,
                              conv_stride=args.stride,
                              hidden_fc_layers_out=args.fcout,
                              fc_dropout=args.dropout,
                              kernel_size=args.kernel)

        lr = args.l
        if args.w:
            weights = torch.FloatTensor(train.get_weights())
            criterion = nn.BCEWithLogitsLoss(pos_weight=weights)
        else:
            criterion = nn.BCEWithLogitsLoss()

        optimizer = Adam(params=model.parameters(), lr=lr)

        scheduler = lr_scheduler.ReduceLROnPlateau(
            optimizer=optimizer,
            mode='min',
            factor=0.5,
            patience=5,
            min_lr=0.00001
        )

        epochs, train_losses, test_losses, accuracies, learning_rates, precisions, recalls, class_precisions, \
            class_recalls, train_accuracies = [], [], [], [], [], [], [], [], [], []

        num_epochs = args.e
        aa_metrics_folder = os.path.join(output_dir, "amino_acid_metrics")
        if not os.path.exists(aa_metrics_folder):
            os.mkdir(aa_metrics_folder)
        aa_metrics = {}
        for epoch in range(num_epochs):
            aa_metrics[epoch + 1] = {}
            for spec in train.all_specificities:
                aa_metrics[epoch + 1][spec] = SubstrateMetrics(spec)

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
                epoch=epoch + 1,
                metrics=aa_metrics,
                specificity_list=train.all_specificities
            )
            scheduler.step(test_loss)
            for aa, metrics in aa_metrics[epoch + 1].items():
                metrics.print()

            epochs.append(epoch + 1)
            train_losses.append(train_loss)
            test_losses.append(test_loss)
            learning_rates.append(lr)
            accuracies.append(accuracy)
            recalls.append(recall)
            class_recalls.append(class_recall)
            precisions.append(precision)
            class_precisions.append(class_precision)
            train_accuracies.append(train_accuracy)

            write_aa_metrics(aa_metrics_folder, epoch + 1, aa_metrics[epoch + 1], train.all_specificities)

            print(
                f"Epoch: {epoch + 1}/{num_epochs}; "
                f"LR: {lr}; "
                f"Train loss: ", train_loss,
                f"Test loss: ", test_loss,
                f"Accuracy: ", accuracy,
                f"Recall: ", recall,
                f"Precision: ", precision,
            )

        out_plots = os.path.join(output_dir, "plots")

        if not os.path.exists(out_plots):
            os.mkdir(out_plots)

        metric_plots = os.path.join(out_plots, "amino_acid_metrics")
        if not os.path.exists(metric_plots):
            os.mkdir(metric_plots)

        plot_aa_metrics(metric_plots, aa_metrics, train.all_specificities)
        plot_metrics(out_plots, epochs, train_losses, test_losses, accuracies, precisions, recalls,
                     class_precisions, class_recalls, train_accuracies)

        metrics_out = os.path.join(output_dir, "metrics.txt")
        write_metrics(metrics_out, epochs, train_losses, test_losses, accuracies, learning_rates, precisions, recalls,
                      class_precisions, class_recalls, train_accuracies)

        save_model(num_epochs - 1, model, optimizer, output_dir)


if __name__ == "__main__":
    crossvalidation()
