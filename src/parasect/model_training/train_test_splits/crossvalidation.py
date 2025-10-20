import os
from collections import Counter
from typing import Optional, Union

from sklearn.model_selection import StratifiedKFold
from iterstrat.ml_stratifiers import MultilabelStratifiedKFold

from parasect.database.build_database import AdenylationDomain, Substrate
from parasect.core.writers import write_list
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode


def make_crossval_sets(train_x: list[AdenylationDomain], train_y: Union[list[Substrate], list[list[Substrate]]],
                       out_dir: str,
                       fold_validation: int = 3, binary_labels: Optional[list] = None,
                       selection_mode: SubstrateSelectionMode = SubstrateSelectionMode.FIRST_ONLY) -> None:
    if selection_mode != SubstrateSelectionMode.ALL:
        stratifier = StratifiedKFold(n_splits=fold_validation, random_state=25051989, shuffle=True)
        splits = stratifier.split(train_x, [s.name for s in train_y])
    else:
        stratifier = MultilabelStratifiedKFold(n_splits=fold_validation, random_state=25051989, shuffle=True)
        splits = stratifier.split(train_x, binary_labels)

    for i, (train_indices, test_indices) in enumerate(splits):
        out_train = os.path.join(out_dir, f"train_crossvalidation_{i + 1}.txt")
        out_test = os.path.join(out_dir, f"test_crossvalidation_{i + 1}.txt")
        counts_out = os.path.join(out_dir, f"crossvalidation_substrate_counts_{i + 1}.txt")

        train = [train_x[j].get_name() for j in train_indices]
        test = [train_x[j].get_name() for j in test_indices]

        if selection_mode != SubstrateSelectionMode.ALL:
            substrates_train = Counter([train_y[j].name for j in train_indices])
            substrates_test = Counter([train_y[j].name for j in test_indices])
        else:
            substrates_train = Counter([substrate.name for j in train_indices for substrate in train_y[j]])
            substrates_test = Counter([substrate.name for j in test_indices for substrate in train_y[j]])

        with open(counts_out, 'w') as out:
            out.write(f"substrate\ttrain\ttest\n")
            for substrate in sorted(set(substrates_train) | set(substrates_test)):
                out.write(f"{substrate}\t{substrates_train[substrate]}\t{substrates_test[substrate]}\n")

        write_list(train, out_train)
        write_list(test, out_test)
