import os


def iterate_over_dir(directory, extension=None, get_dirs=False):
    for file_name in os.listdir(directory):
        if not get_dirs:
            if file_name.endswith(extension):
                file_label = file_name.split(extension)[0]
                file_path = os.path.join(directory, file_name)
                yield file_label, file_path
        else:
            file_path = os.path.join(directory, file_name)
            if os.path.isdir(file_path):
                yield file_name, file_path


def find_crossval_pairs(input_dir):
    crossval_set_to_data = {}
    for input_file in os.listdir(input_dir):
        input_path = os.path.join(input_dir, input_file)
        if os.path.isfile(input_path) and input_file.endswith('.txt'):
            _, data_type, crossval_set = input_file[:-4].split('_')
            crossval_set = int(crossval_set)
            if crossval_set not in crossval_set_to_data:
                crossval_set_to_data[crossval_set] = {"train": '',
                                                      "test": ''}
            crossval_set_to_data[crossval_set][data_type] = input_path

    path_pairs = []

    for crossval_set, type_to_path in crossval_set_to_data.items():
        assert type_to_path["train"] and type_to_path["test"]
        path_pairs.append((type_to_path["train"], type_to_path["test"]))

    path_pairs.sort()
    return path_pairs