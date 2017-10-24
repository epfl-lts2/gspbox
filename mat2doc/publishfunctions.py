import os
import commands


def dos2unix(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name = os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('dos2unix '+name)


def unix2dos(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name = os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('unix2dos '+name)
