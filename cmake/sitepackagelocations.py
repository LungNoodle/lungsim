#!/usr/bin/env python

import os, sys, site

PREFIXES = [sys.prefix, sys.exec_prefix]


def makepath(*paths):
    dir = os.path.join(*paths)
    try:
        dir = os.path.abspath(dir)
    except OSError:
        pass
    return dir, os.path.normcase(dir)


def local_getsitepackages():
    sitepackages = []
    seen = set()

    for prefix in PREFIXES:
        if not prefix or prefix in seen:
            continue
        seen.add(prefix)

        if sys.platform in ('os2emx', 'riscos'):
            sitepackages.append(os.path.join(prefix, "Lib", "site-packages"))
        elif os.sep == '/':
            sitepackages.append(os.path.join(prefix, "lib",
                                        "python" + sys.version[:3],
                                        "site-packages"))
            sitepackages.append(os.path.join(prefix, "lib", "site-python"))
        else:
            sitepackages.append(prefix)
            sitepackages.append(os.path.join(prefix, "lib", "site-packages"))
        if sys.platform == "darwin":
            # for framework builds *only* we add the standard Apple
            # locations.
            from sysconfig import get_config_var
            framework = get_config_var("PYTHONFRAMEWORK")
            if framework:
                sitepackages.append(
                        os.path.join("/Library", framework,
                            sys.version[:3], "site-packages"))
    return sitepackages


if hasattr(sys, 'real_prefix'):
    # This is a virtual environment and probably doesn't have getsitepackages() for site.
    if not hasattr(site, 'getsitepackages'):
        site.getsitepackages = local_getsitepackages

site_packages = site.getsitepackages()
if sys.platform == "darwin":
    known_paths = []
    current_site_packages = site_packages[:]
    for p in current_site_packages:
        fullname = os.path.join(p, 'Extras.pth')
        try:
            f = open(fullname, "rU")
        except IOError:
            continue
        with f:
            for n, line in enumerate(f):
                if line.startswith("#"):
                    continue
                try:
                    if line.startswith(("import ", "import\t")):
                        continue
                    line = line.rstrip()
                    dir, dircase = makepath(line, '')
                    if not dircase in known_paths and os.path.exists(dir):
                        site_packages.append(dir)
                        known_paths.add(dircase)
                except Exception:
                    pass


for site_package in site_packages:
    sys.stdout.write('%s\n' % site_package)
