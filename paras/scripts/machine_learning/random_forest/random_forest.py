import os

from joblib import dump, load
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE
from matplotlib import pyplot as plt
from numpy import argmax, sqrt

from paras.scripts.data_analysis.plotting.plot_confusion_matrices import plot_matrix, write_matrix
from statistics import mean, median, stdev


class DataPoint:
    def __init__(self, domain: str):
        self.domain = domain
        self.domain_substrates = []
        self.substrate = None
        self.vector = []
        self.y = None

    def extend_vector(self, vector_mapping, data_type='domain'):
        if data_type == 'substrate':
            self.vector += vector_mapping[self.substrate]
        elif data_type == 'domain':
            self.vector += vector_mapping[self.domain]
        else:
            raise ValueError("Data type must be 'substrate' or 'domain'")

    def set_y(self, y):
        self.y = y


def get_mean(numbers):
    if numbers:
        return mean(numbers)
    else:
        return "N/A"


def get_stdev(numbers):
    if numbers and len(numbers) > 1:
        return stdev(numbers)
    else:
        return "N/A"


def get_median(numbers):
    if numbers:
        return median(numbers)
    else:
        return "N/A"


class DomainDataPoint(DataPoint):
    def __init__(self, domain):
        super().__init__(domain)


class TandemDataPoint(DataPoint):
    def __init__(self, domain, compound):
        super().__init__(domain)
        self.substrate = compound


class Dataset:
    def __init__(self, domain_list, train_domain_list, test_domain_list, substrate_list, domain_to_substrate,
                 domain_mappings, domain_mapping_categories, compound_mappings=None,
                 compound_mapping_categories=None, mode='single_label'):

        self.train_datapoints = []
        self.test_datapoints = []
        self.substrates = substrate_list
        self.categories = []
        self.domain_to_substrate = domain_to_substrate

        for mapping_type in domain_mappings:
            self.categories += domain_mapping_categories[mapping_type]

        self.type = mode

        if self.type == 'tandem':
            assert compound_mappings is not None
            assert compound_mapping_categories is not None

            for mapping_type in compound_mappings:
                self.categories += compound_mapping_categories[mapping_type]

        for domain in domain_list:

            if domain in train_domain_list or domain in test_domain_list:
                specs = domain_to_substrate[domain]

                substrates = []

                for substrate in specs:
                    if substrate in self.substrates:
                        substrates.append(substrate)

                if compound_mappings:
                    for substrate in self.substrates:
                        datapoint = TandemDataPoint(domain, substrate)
                        for mapping_type, mapping in domain_mappings.items():
                            datapoint.extend_vector(mapping)
                        for mapping_type, mapping in compound_mappings.items():
                            datapoint.extend_vector(mapping, data_type='substrate')

                        if substrate in substrates:
                            datapoint.y = 1
                        else:
                            datapoint.y = 0

                        datapoint.domain_substrates = substrates

                        if domain in test_domain_list:
                            self.test_datapoints.append(datapoint)
                        elif domain in train_domain_list:
                            self.train_datapoints.append(datapoint)

                else:
                    datapoint = DomainDataPoint(domain)
                    datapoint.domain_substrates = substrates
                    for mapping_type, mapping in domain_mappings.items():
                        datapoint.extend_vector(mapping)

                    if mode == 'single_label':
                        datapoint.y = substrates[0]
                    elif mode == 'multi_label':
                        datapoint.y = substrates
                    elif mode == 'binary':
                        datapoint.y = []
                        for substrate in self.substrates:
                            if substrate in substrates:
                                datapoint.y.append(1)
                            else:
                                datapoint.y.append(0)
                    else:
                        raise ValueError("Only 'single_label', 'multi_label' and 'binary' are accepted modes for domain-only input.")

                    if domain in test_domain_list:
                        self.test_datapoints.append(datapoint)
                    elif domain in train_domain_list:
                        self.train_datapoints.append(datapoint)

    def get_data(self, mode):
        if mode == 'train':
            return self.train_datapoints
        elif mode == 'test':
            return self.test_datapoints
        else:
            raise ValueError("Data retrieval mode must be 'train' or 'test'.")

    def get_vectors(self, mode='train'):
        if mode == 'train':
            datapoints = self.train_datapoints
        elif mode == 'test':
            datapoints = self.test_datapoints
        else:
            raise ValueError("Vector retrieval mode must be 'train' or 'test'.")

        x = []
        y = []

        for datapoint in datapoints:
            x.append(datapoint.vector)
            y.append(datapoint.y)

        return x, y

    def write_test_metrics(self, mode, predictions, out_dir, oob_score=None):
        if mode == 'train':
            datapoints = self.train_datapoints
            assert oob_score is not None
        elif mode == 'test':
            datapoints = self.test_datapoints
        else:
            raise ValueError("Data retrieval mode must be 'train' or 'test'.")

        confusion_matrix_file = os.path.join(out_dir, "confusion_matrix.txt")
        confusion_matrix_plot = os.path.join(out_dir, "confusion_matrix.svg")
        metrics_file = os.path.join(out_dir, "metrics.txt")
        accuracy_file = os.path.join(out_dir, "accuracy.txt")
        test_results_file = os.path.join(out_dir, "test_results.txt")
        interaction_probabilities = os.path.join(out_dir, "interaction_probabilities.txt")

        correct = 0
        incorrect = 0

        correct_probabilities = []
        incorrect_probabilities = []

        substrate_to_metrics = {}

        confusion_matrix = []

        tp_probabilities = []
        fp_probabilities = []
        fn_probabilities = []
        tn_probabilities = []

        fp = 0
        fn = 0
        tp = 0
        tn = 0

        for substrate in self.substrates:
            row = []
            substrate_to_metrics[substrate] = {}
            substrate_to_metrics[substrate]["correct"] = 0
            substrate_to_metrics[substrate]["incorrect"] = 0

            substrate_to_metrics[substrate]["TP"] = 0
            substrate_to_metrics[substrate]["FP"] = 0
            substrate_to_metrics[substrate]["FN"] = 0
            substrate_to_metrics[substrate]["TN"] = 0

            substrate_to_metrics[substrate]["TP probabilities"] = []
            substrate_to_metrics[substrate]["FP probabilities"] = []
            substrate_to_metrics[substrate]["FN probabilities"] = []
            substrate_to_metrics[substrate]["TN probabilities"] = []
            substrate_to_metrics[substrate]["correct probabilities"] = []
            substrate_to_metrics[substrate]["incorrect probabilities"] = []

            for _ in self.substrates:
                row.append(0)

            confusion_matrix.append(row)

        domain_to_top_prediction = {}

        if self.type == 'single_label':

            for i, datapoint in enumerate(datapoints):
                probability, prediction = predictions[i][0]
                domain_to_top_prediction[datapoint.domain] = (prediction, probability)

                if prediction in datapoint.domain_substrates:
                    correct += 1
                    substrate_to_metrics[datapoint.y]["correct"] += 1
                    substrate_to_metrics[datapoint.y]["TP"] += 1
                    substrate_to_metrics[datapoint.y]["TP probabilities"].append(probability)
                    tp_probabilities.append(probability)
                else:
                    incorrect += 1
                    substrate_to_metrics[datapoint.y]["incorrect"] += 1
                    substrate_to_metrics[datapoint.y]["FN"] += 1
                    for prob, pred in predictions[i]:
                        if pred == datapoint.y:
                            substrate_to_metrics[datapoint.y]["FN probabilities"].append(prob)
                            fn_probabilities.append(prob)

                    substrate_to_metrics[prediction]["FP"] += 1
                    substrate_to_metrics[prediction]["FP probabilities"].append(probability)
                    fp_probabilities.append(probability)

                confusion_matrix[self.substrates.index(prediction)][self.substrates.index(datapoint.y)] += 1

        elif self.type == 'tandem':

            domain_to_predictions = {}
            domain_to_correct = {}
            domain_to_best_prob = {}

            for i, datapoint in enumerate(datapoints):
                true_y = datapoint.y
                prob_per_prediction = predictions[i]
                probability_predicted_y, predicted_y = prob_per_prediction[0]

                if predicted_y:
                    prob_true = prob_per_prediction[0][0]
                else:
                    prob_true = prob_per_prediction[1][0]

                if datapoint.domain not in domain_to_correct:
                    domain_to_correct[datapoint.domain] = False
                    domain_to_predictions[datapoint.domain] = []
                    domain_to_best_prob[datapoint.domain] = 0.0

                if prob_true > domain_to_best_prob[datapoint.domain]:
                    domain_to_best_prob[datapoint.domain] = prob_true
                    if datapoint.substrate in datapoint.domain_substrates:
                        domain_to_correct[datapoint.domain] = True
                    else:
                        domain_to_correct[datapoint.domain] = False

                domain_to_predictions[datapoint.domain].append((prob_true, datapoint.substrate,
                                                                datapoint.substrate in datapoint.domain_substrates,
                                                                datapoint))

                if true_y:
                    if not predicted_y:
                        fn += 1
                        fn_probabilities.append(prob_true)
                        substrate_to_metrics[datapoint.substrate]["FN probabilities"].append(prob_true)
                    else:
                        tp += 1
                        tp_probabilities.append(prob_true)
                        substrate_to_metrics[datapoint.substrate]["TP probabilities"].append(prob_true)
                        confusion_matrix[self.substrates.index(datapoint.substrate)][
                            self.substrates.index(datapoint.substrate)] += 1

                elif not true_y:
                    if not predicted_y:
                        tn += 1
                        tn_probabilities.append(prob_true)
                        substrate_to_metrics[datapoint.substrate]["TN probabilities"].append(prob_true)
                    else:
                        fp += 1
                        fp_probabilities.append(prob_true)
                        substrate_to_metrics[datapoint.substrate]["FP probabilities"].append(prob_true)

                        for true_substrate in datapoint.domain_substrates:
                            confusion_matrix[self.substrates.index(datapoint.substrate)][self.substrates.index(true_substrate)] += 1

            correct = list(domain_to_correct.values()).count(True)
            incorrect = list(domain_to_correct.values()).count(False)

            with open(interaction_probabilities, 'w') as probabilities_file:
                probabilities_file.write("domain\tsubstrate\tinteraction_probability\ttrue_interaction\n")

                for domain in domain_to_correct:
                    predictions = domain_to_predictions[domain]

                    predictions.sort(reverse=True)

                    for prediction in predictions:
                        probability = prediction[0]
                        substrate = prediction[1]
                        true_interaction = prediction[2]
                        datapoint = prediction[3]
                        probabilities_file.write(
                            f"{datapoint.domain}\t{substrate}\t{probability:.2f}\t{true_interaction}\n")

                    probability = predictions[0][0]
                    datapoint = predictions[0][3]
                    domain_to_top_prediction[domain] = (predictions[0][1], probability)

                    for substrate in datapoint.domain_substrates:
                        if domain_to_correct[domain]:
                            substrate_to_metrics[substrate]["correct probabilities"].append(probability)
                            substrate_to_metrics[substrate]["correct"] += 1
                            correct_probabilities.append(probability)
                        else:
                            substrate_to_metrics[substrate]["incorrect probabilities"].append(probability)
                            substrate_to_metrics[substrate]["incorrect"] += 1
                            incorrect_probabilities.append(probability)

        with open(metrics_file, 'w') as out:
            if self.type != 'tandem':
                out.write("substrate\taccuracy\tnr_correct\tnr_incorrect\tprecision\trecall\tF1\tTP\tTP_probability\tFP\tFP_probability\tFN\tFN_probability\n")
            else:
                out.write(
                    "substrate\taccuracy\tnr_correct\tcorrect_probability\tnr_incorrect\tincorrect_probability\tprecision\trecall\tF1\tTP\tTP_probability\tFP\tFP_probability\tFN\tFN_probability\tTN\tTN_probability\n")

            for substrate in self.substrates:
                nr_tp = substrate_to_metrics[substrate]["TP"]
                nr_fp = substrate_to_metrics[substrate]["FP"]
                nr_fn = substrate_to_metrics[substrate]["FN"]
                tp_prob = get_mean(substrate_to_metrics[substrate]["TP probabilities"])
                fp_prob = get_mean(substrate_to_metrics[substrate]["FP probabilities"])
                fn_prob = get_mean(substrate_to_metrics[substrate]["FN probabilities"])

                nr_correct = substrate_to_metrics[substrate]["correct"]
                nr_incorrect = substrate_to_metrics[substrate]["incorrect"]
                if nr_correct + nr_incorrect == 0:
                    accuracy = "N/A"
                else:
                    accuracy = nr_correct / (nr_incorrect + nr_correct)

                if nr_tp + nr_fp == 0:
                    precision = 0.0
                else:
                    precision = nr_tp / (nr_tp + nr_fp)
                if nr_tp + nr_fn == 0:
                    recall = 0.0
                else:

                    recall = nr_tp / (nr_tp + nr_fn)
                if precision + recall == 0:
                    f1 = 0.0
                else:
                    f1 = 2 * precision * recall / (precision + recall)
                if self.type != 'tandem':

                    out.write(
                        f"{substrate}\t{accuracy}\t{nr_correct}\t{nr_incorrect}\t{precision}\t{recall}\t{f1}\t{nr_tp}\t{tp_prob}\t{nr_fp}\t{fp_prob}\t{nr_fn}\t{fn_prob}\n")
                else:
                    nr_tn = substrate_to_metrics[substrate]["TN"]
                    tn_prob = get_mean(substrate_to_metrics[substrate]["TN probabilities"])

                    correct_prob = get_mean(substrate_to_metrics[substrate]["correct probabilities"])
                    incorrect_prob = get_mean(substrate_to_metrics[substrate]["incorrect probabilities"])
                    out.write(
                        f"{substrate}\t{accuracy}\t{nr_correct}\t{correct_prob}\t{nr_incorrect}\t{incorrect_prob}\t{precision}\t{recall}\t{f1}\t{nr_tp}\t{tp_prob}\t{nr_fp}\t{fp_prob}\t{nr_fn}\t{fn_prob}\t{nr_tn}\t{tn_prob}\n")
        box_out = os.path.join(out_dir, "distributions.svg")
        if self.type == 'single_label':

            boxplot = plt.boxplot([tp_probabilities, fp_probabilities, fn_probabilities],
                                  vert=True,
                                  labels=["TP", "FP", "FN"],
                                  patch_artist=True)
            colours = ["lightsteelblue", "navajowhite", "pink"]
            ylabel = "Prediction confidence"
        elif self.type == 'tandem':
            boxplot = plt.boxplot([tp_probabilities, fp_probabilities, fn_probabilities, tn_probabilities],
                                  vert=True,
                                  labels=["TP", "FP", "FN", "TN"],
                                  patch_artist=True)
            colours = ["lightsteelblue", "navajowhite", "pink", "lightsteelblue"]
            ylabel = "Predicted interaction probability"
        else:
            raise ValueError("'tandem' and 'single_label' are the only supported modes.")
        for i, patch in enumerate(boxplot['boxes']):
            colour = colours[i]
            patch.set_facecolor(colour)
        plt.grid(True)
        plt.ylabel(ylabel)
        plt.savefig(box_out, bbox_inches='tight')

        with open(accuracy_file, 'w') as out:
            accuracy = correct / (correct + incorrect)
            mean_tp = get_mean(tp_probabilities)
            mean_fp = get_mean(fp_probabilities)
            mean_fn = get_mean(fn_probabilities)

            median_tp = get_median(tp_probabilities)
            median_fp = get_median(fp_probabilities)
            median_fn = get_median(fn_probabilities)

            stdev_tp = get_stdev(tp_probabilities)
            stdev_fp = get_stdev(fp_probabilities)
            stdev_fn = get_stdev(fn_probabilities)

            out.write(f"accuracy\t{accuracy}\n")
            out.write(f"correct\t{correct}\n")
            out.write(f"incorrect\t{incorrect}\n")
            out.write(f"mean_probability_tp\t{mean_tp}\n")
            out.write(f"mean_probability_fp\t{mean_fp}\n")
            out.write(f"mean_probability_fn\t{mean_fn}\n")
            out.write(f"median_probability_tp\t{median_tp}\n")
            out.write(f"median_probability_fp\t{median_fp}\n")
            out.write(f"median_probability_fn\t{median_fn}\n")
            out.write(f"stdev_probability_tp\t{stdev_tp}\n")
            out.write(f"stdev_probability_fp\t{stdev_fp}\n")
            out.write(f"stdev_probability_fn\t{stdev_fn}\n")

            if self.type == 'tandem':
                mean_correct = get_mean(correct_probabilities)
                mean_incorrect = get_mean(incorrect_probabilities)
                median_correct = get_median(correct_probabilities)
                median_incorrect = get_median(incorrect_probabilities)
                stdev_correct = get_stdev(correct_probabilities)
                stdev_incorrect = get_stdev(incorrect_probabilities)

                out.write(f"mean_probability_correct\t{mean_correct}\n")
                out.write(f"median_probability_correct\t{median_correct}\n")
                out.write(f"mean_probability_incorrect\t{mean_incorrect}\n")
                out.write(f"median_probability_incorrect\t{median_incorrect}\n")
                out.write(f"stdev_probability_correct\t{stdev_correct}\n")
                out.write(f"stdev_probability_incorrect\t{stdev_incorrect}\n")

                out.write(f"true positives\t{tp}\n")
                out.write(f"false positives\t{fp}\n")
                out.write(f"true negatives\t{tn}\n")
                out.write(f"false negatives\t{fn}\n")

            if oob_score:
                out.write(f"oob_score\t{oob_score}\n")

        with open(test_results_file, 'w') as out:
            out.write("domain\tprediction\tconfidence\n")
            for domain, top_prediction in domain_to_top_prediction.items():
                prediction, probability = top_prediction
                out.write(f"{domain}\t{prediction}\t{probability}\n")

        plot_matrix(confusion_matrix, self.substrates, confusion_matrix_plot)
        write_matrix(confusion_matrix, self.substrates, confusion_matrix_file)


class RandomForest:
    def __init__(self, sampling_method=None, nr_trees=1000, threads=-1, model=None):
        if model is None:
            self.sampling_method = sampling_method
            if self.sampling_method == 'balanced':

                classifier = RandomForestClassifier(n_estimators=nr_trees, n_jobs=threads, oob_score=True,
                                                    random_state=25051989, class_weight='balanced')
            elif self.sampling_method == 'balanced_subsample':
                classifier = RandomForestClassifier(n_estimators=nr_trees, n_jobs=threads, oob_score=True,
                                                    random_state=25051989, class_weight='balanced_subsample')
            else:
                classifier = RandomForestClassifier(n_estimators=nr_trees, n_jobs=threads, oob_score=True,
                                                    random_state=25051989)

            self.classifier = classifier
        else:
            self.sampling_method = sampling_method
            self.classifier = load(model)

    def train(self, dataset, out_path=None):
        train_x, train_y = self.build_data(dataset, mode='train', apply_sampling=True)
        print("Training classifier..")
        self.classifier.fit(train_x, train_y)
        print(f"Out of bag score: {self.classifier.oob_score_}")
        if out_path:
            dump(self.classifier, out_path)

    def test(self, dataset, mode='test'):
        test_x, test_y = self.build_data(dataset, mode)
        classes = self.classifier.classes_
        predict_y = self.classifier.predict_proba(test_x)

        predictions = []
        for prediction in predict_y:
            prob_per_class = list(zip(prediction, classes))

            prob_per_class.sort(reverse=True)
            predictions.append(prob_per_class)

        return predictions, test_y, predict_y

    def record_test_results(self, dataset, out_dir, mode='both'):
        if mode == 'both' or mode == 'train':
            train_predictions, test_y, predict_y = self.test(dataset, mode='train')
            train_performance = os.path.join(out_dir, 'train_performance')
            if not os.path.exists(train_performance):
                os.mkdir(train_performance)
            dataset.write_test_metrics("train", train_predictions, train_performance, self.classifier.oob_score_)
        if mode == 'both' or mode == 'test':
            test_predictions, test_y, predict_y = self.test(dataset, mode='test')

            test_performance = os.path.join(out_dir, 'test_performance')
            if not os.path.exists(test_performance):
                os.mkdir(test_performance)

            if dataset.type == 'tandem':

                precision, recall, thresholds = precision_recall_curve(test_y, predict_y[:, 1])
                fscore = (2 * precision * recall) / (precision + recall)
                index = argmax(fscore)

                plt.plot(recall, precision)
                plt.scatter(recall[index], precision[index], marker='o', color='black', label='Recommended threshold')
                plt.text(recall[index] - 0.65, precision[index] - 0.05, f"Recommended threshold: {thresholds[index]:.2f}",
                         size=12,
                         fontdict=dict(size=10),
                         bbox=dict(facecolor="lightsteelblue", alpha=0.5)
                         )
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                plt.savefig(os.path.join(test_performance, "pr_curve.svg"), bbox_inches='tight')
                plt.clf()

                fpr, tpr, threshold = roc_curve(test_y, predict_y[:, 1])

                gmeans = sqrt(tpr * (1 - fpr))
                ix = argmax(gmeans)
                plt.plot(fpr, tpr)
                plt.scatter(fpr[ix], tpr[ix], marker='o', color='black', label='Recommended threshold')
                plt.text(fpr[ix] + 0.05, tpr[ix] - 0.05, f"Recommended threshold: {thresholds[ix]:.2f}",
                         size=12,
                         fontdict=dict(size=10),
                         bbox=dict(facecolor="lightsteelblue", alpha=0.5)
                         )
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.savefig(os.path.join(test_performance, "roc_curve.svg"), bbox_inches='tight')
                plt.clf()

            dataset.write_test_metrics("test", test_predictions, test_performance)

    def build_data(self, dataset, mode, apply_sampling=False):
        x, y = dataset.get_vectors(mode)

        if apply_sampling:

            if self.sampling_method == 'under_sample':
                print("Under-sampling..")

                undersampler = RandomUnderSampler(random_state=25051989)

                under_x, under_y = undersampler.fit_resample(x, y)

                return under_x, under_y

            elif self.sampling_method == 'over_sample':
                print("Over-sampling..")
                oversampler = SMOTE(random_state=25051989)
                over_x, over_y = oversampler.fit_resample(x, y)
                return over_x, over_y

            else:
                return x, y

        else:
            return x, y

    def write_feature_importances(self, categories, out_dir):
        out_file = os.path.join(out_dir, 'feature_importances.txt')
        assert len(categories) == len(self.classifier.feature_importances_)
        importance_and_categories = []
        with open(out_file, 'w') as out:
            out.write("category\timportance\n")

            for i, feature_importance in enumerate(self.classifier.feature_importances_):
                category = categories[i]
                importance_and_categories.append((feature_importance, category))

            importance_and_categories.sort(reverse=True)
            for importance, category in importance_and_categories:
                out.write(f"{category}\t{importance:.5f}\n")

