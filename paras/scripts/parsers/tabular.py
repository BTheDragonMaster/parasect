from collections import OrderedDict


class Tabular:
    def __init__(self, tabular_path, id_columns, separator='\t'):
        self.separator = separator
        self.data = OrderedDict()
        with open(tabular_path, 'r') as tabular_file:
            self.header = tabular_file.readline()
            header = self.header.split(self.separator)
            self.categories = []
            for i, category in enumerate(header):
                self.categories.append(category.strip())
            for line in tabular_file:

                values = line.split(self.separator)
                row_id = []
                for id_column in id_columns:
                    row_id.append(values[id_column])

                row_id = tuple(row_id)

                if row_id in self.data:
                    print(f"WARNING! Duplicate row ID found: {row_id}")
                self.data[row_id] = OrderedDict()

                for j, value in enumerate(values):
                    category = self.categories[j]
                    self.data[row_id][category] = value.strip()

    def get_column(self, category):
        if category not in self.categories:
            raise KeyError(f"Cannot find category {category} in data.")

        column = []

        for data_id in self.data:
            column.append(self.get_value(data_id, category))

        return column

    def get_row(self, data_id):
        if data_id not in self.data:
            raise KeyError(f"Cannot find data ID {data_id} in data.")

        row = []

        for category in self.categories:
            row.append(self.get_value(data_id, category))

        return row

    def get_value(self, data_id, category):
        try:
            return self.data[data_id][category]
        except KeyError:
            if data_id in self.data:
                print(f"Cannot find category {category} in data.")
            else:
                print(f"Cannot find id {data_id} in data.")
            raise KeyError

    def write_table(self, out_file):
        with open(out_file, 'w') as out:
            out.write(self.header)
            for seq_id in self.data:
                for i, category in enumerate(self.categories):
                    if i == len(self.categories) - 1:
                        out.write(f"{self.data[seq_id][category]}\n")
                    else:
                        out.write(f"{self.data[seq_id][category]}{self.separator}")
