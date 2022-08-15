
=========
Acquiring
=========

The first step to working with the pulmonary simulation library is to acquire it.
The easiest way to obtain the codebase is to use git and clone the source code.
It is also possible to acquire the codebase as a zip archive if git is not available.

---
Git
---

Git is the distributed source control management tool in wide use by millions of software developers and engineers.
If git is not available on your platform the following sub-sections provide guidance an how to acquire it.

Windows
=======

On Windows download `Git for windows <https://git-scm.com/download/win>`_ and install in the usual manner.

OS X
====

OS X comes preloaded with git and is available directly from the command prompt.

GNU/Linux
=========

The package manager for the distro will (most likely) have the required package to install.
Before installing check to see if git is already available::

  git --version
  
If git is not available you can install the package on an Ubuntu distribution with the following command::

  sudo apt-get install git

Other distributions will have a similar command.

Get the Source Code with Git
============================

With git installed we can acquire the codebase with the command::

  git clone https://github.com/LungNoodle/lungsim.git

---
Zip
---

To download a zip archive of the current *develop* branch use a web browser to visit this page https://github.com/LungNoodle/lungsim/archive/develop.zip.
Then simply unzip the archive::

  unzip lungsim-develop.zip
