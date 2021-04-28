

import os


def imgpath(fn):
    if not os.path.exists('img'):
        os.mkdir('img')
    return '/'.join(['img', fn])


def showfigs():
    return True


