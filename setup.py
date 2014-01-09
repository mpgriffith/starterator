from ez_setup import use_setuptools
use_setuptools()

import os, shutil, sys, subprocess, ConfigParser
import setuptools
def desktop_file():
    # if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "starterator.svg")):
    shutil.copyfile("starterator/extras/starterator.svg",
            os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.svg"))
    # if not os.path.exists(os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.desktop")):
    config = ConfigParser.RawConfigParser()
    desktop = config.read("starterator/extras/starterator.desktop")
    # desktop.set("")
    shutil.copyfile("starterator/extras/starterator.desktop",
            os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.desktop"))

def install_pyPDF():
    print "installing PyPDF2"
    subprocess.check_call(["pip", "install", "PyPDF2"])
    print "PyPDF2 installed"
    try:
        import PyPDF2
        print "yay"
    except ImportError:
        print "boo"

if "install" in sys.argv:
    desktop_file()
    # install_pyPDF()

setuptools.setup(  
    name='Starterator',
    version='1.0',
    author='Marissa Pacey',
    author_email='mppacey@gmail.com',
    packages=setuptools.find_packages(),
    install_requires=["PyPDF2", "MySQl-python", "reportlab"],
    entry_points={"gui_scripts": ["Starterator = starterator.window:gui"]},
    # package_dir = {"":"src"},
    package_data={"": ["extras/*", "Help/*.page", "Help/figures/*"]},
    # data_files={os.environ["HOME"]+"/.local/share/applications":
    #                 ["starterator/extras/starterator.desktop", "staterator/extras/starterator_icon.svg"],
    #             os.environ["HOME"]+"/Starterator/": ["Help/*", "starterator/extras/starterator.config"]
    #  },
    description='Starterator - a program to compare start sites in phamilies'
)

