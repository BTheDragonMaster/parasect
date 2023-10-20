import os

import paras.tmp
import shutil

TEMP_DIR = os.path.dirname(paras.tmp.__file__)


def clear_temp():

    for temp_name in os.listdir(TEMP_DIR):
        if not temp_name == '__init__.py':
            file_location = os.path.join(TEMP_DIR, temp_name)
            if os.path.isfile(file_location):
                os.remove(file_location)
            elif os.path.isdir(file_location):
                shutil.rmtree(file_location)
