configure_file(@CMAKE_CURRENT_SOURCE_DIR@/__init__.py ${PYTHON_PACKAGE_DIR}/__init__.py)
configure_file(@AETHER_README_FILE@ ${PYTHON_PACKAGE_DIR}/../README.rst)
configure_file(@SETUP_PRE_GEN_PY_FILE@ ${PYTHON_PACKAGE_DIR}/../setup.py)
