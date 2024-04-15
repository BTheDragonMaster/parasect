import os
import shutil

def clear_temp(temp_dir):

    for temp_name in os.listdir(temp_dir):
        
        # skip .gitkeep 
        if not temp_name == '__init__.py' and not temp_name == '.gitkeep':
            file_location = os.path.join(temp_dir, temp_name)
            if os.path.isfile(file_location):
                os.remove(file_location)
            elif os.path.isdir(file_location):
                shutil.rmtree(file_location)
