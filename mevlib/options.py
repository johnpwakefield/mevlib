

#TODO get rid of this file


import os


def imgpath(fn):
    if not os.path.exists('img'):
        os.mkdir('img')
    return '/'.join(['img', fn])


def showfigs():
    return False


