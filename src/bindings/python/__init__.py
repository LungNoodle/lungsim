import os
import platform
import sys


if platform.system() == "Windows":
    os.add_dll_directory(os.path.join(sys.prefix, "Library", "bin"))

