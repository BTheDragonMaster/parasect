from sys import argv
import os
from typing import Any

import numpy as np

from parasect.core.tabular import Tabular, write_tabular
from parasect.core.parsing import iterate_over_dir


from pathlib import Path

def find_best_row_file(tabulars: list[Tabular], out_file: str, metric_column: str) -> str:
    """
    Computes the average table, writes it to out_file,
    and returns the row ID corresponding to the best value in metric_column.
    """
    avg_dicts, header = average_tabulars(tabulars)
    write_tabular(avg_dicts, header, out_file)

    # Find best-performing row
    col_index = header.index(metric_column) - 1  # subtract 1 because avg_dicts skips ID column
    column_dict = avg_dicts[col_index]

    best_row_id = max(
        (row_id for row_id, val in column_dict.items() if val != "N/A"),
        key=lambda rid: float(column_dict[rid])
    )
    return best_row_id  # or Path(out_file) if you want the path

def average_tabulars(tabular_list: list[Tabular]) -> tuple[list[dict[str, Any]], list[str]]:
    """
    Compute element-wise averages across multiple Tabular objects, ignoring 'N/A'.

    Returns:
        - A list of dictionaries, one per non-ID column.
          Each dictionary maps row_id -> averaged value.
        - The list of column names from the original Tabular (including ID column).

    Example:
        averaged_dicts, header = average_tabulars([tab1, tab2, tab3])
    """
    if not tabular_list:
        raise ValueError("No Tabular objects provided")

    # Assume all Tabulars have the same structure
    row_ids = list(tabular_list[0].rows.keys())
    column_names = tabular_list[0].column_names

    # Create one dictionary per non-ID column to hold averaged values
    averaged_column_dicts: list[dict[str, Any]] = [dict() for _ in column_names[1:]]

    for row_id in row_ids:
        for col_index, column_name in enumerate(column_names[1:], start=1):  # skip ID column
            numeric_values: list[float] = []

            for tabular in tabular_list:
                cell_value = tabular.get_row_value(row_id, column_name)
                if cell_value != "N/A":
                    try:
                        numeric_values.append(float(cell_value))
                    except ValueError:
                        # Ignore non-convertible entries
                        continue

            if numeric_values:
                average_value = np.mean(numeric_values)
                averaged_column_dicts[col_index - 1][row_id] = str(average_value)
            else:
                averaged_column_dicts[col_index - 1][row_id] = "N/A"

    return averaged_column_dicts, column_names

def main():
    for crossval_category, crossval_category_dir in iterate_over_dir(argv[1], get_dirs=True):

        category_to_parasect = {}
        category_to_paras_confidence = {}
        category_to_paras_accuracy = {}
        for crossval_name, crossval_dir in iterate_over_dir(crossval_category_dir, get_dirs=True):
            category = crossval_name.split('_')[-1]

            if 'first_valid' in crossval_category or 'first_only' in crossval_category:
                if category not in category_to_paras_confidence:
                    category_to_paras_confidence[category] = []
                    category_to_paras_accuracy[category] = []
                confidence = os.path.join(crossval_dir, "confidence.txt")
                accuracy = os.path.join(crossval_dir, "accuracy.txt")
                category_to_paras_confidence[category].append(Tabular(confidence))
                category_to_paras_accuracy[category].append(Tabular(accuracy))
            else:
                if category not in category_to_parasect:
                    parasect_metrics = os.path.join(crossval_dir, "parasect_performance.txt")
                    if category not in category_to_parasect:
                        category_to_parasect[category] = []
                    category_to_parasect[category].append(Tabular(parasect_metrics))
        if category_to_parasect:
            for category, tabulars in category_to_parasect.items():
                average_tabular_dicts, header = average_tabulars(tabulars)
                write_tabular(average_tabular_dicts, header, os.path.join(crossval_category_dir, f"parasect_performance_average_{category}.txt"))
        elif category_to_paras_confidence:
            for category, tabulars in category_to_paras_confidence.items():
                average_tabular_dicts, header = average_tabulars(tabulars)
                write_tabular(average_tabular_dicts, header,
                              os.path.join(crossval_category_dir, f"average_paras_confidence_{category}.txt"))
            for category, tabulars in category_to_paras_accuracy.items():
                average_tabular_dicts, header = average_tabulars(tabulars)
                write_tabular(average_tabular_dicts, header,
                              os.path.join(crossval_category_dir, f"average_paras_accuracy_{category}.txt"))

if __name__ == "__main__":
    main()
