import os
import paras.data.temp

import paras.data.temp

TEMP_DIR = os.path.dirname(paras.data.temp.__file__)


def clear_temp():

    for temp_name in os.listdir(TEMP_DIR):
        if not temp_name == '__init__.py':
            file_location = os.path.join(TEMP_DIR, temp_name)
            if os.path.isfile(file_location):
                os.remove(file_location)