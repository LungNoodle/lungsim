import os
import platform
import sys


if platform.system() == "Windows":
    intel_library_path = os.path.join(sys.prefix, "Library", "bin")
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        os.environ["PATH"] = intel_library_path + os.path.sep + os.environ["PATH"]
    else:
        os.add_dll_directory(intel_library_path)
