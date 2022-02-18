import os
import platform
import sys


if platform.system() == "Windows":
    version_info = sys.version_info
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        sys.path.append(os.path.join(sys.prefix, "Library", "bin"))
    else:
        os.add_dll_directory(os.path.join(sys.prefix, "Library", "bin"))
