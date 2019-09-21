import numpy as np
import datetime
from PIL import Image, ImageChops

def trim(im):
    # trim whitespace from the image file
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)

def yn_choice(message, default='y'):
    choices = 'Y/n' if default.lower() in ('y', 'yes') else 'y/N'
    choice = input("%s (%s) " % (message, choices))
    values = ('y', 'yes', '') if default == 'y' else ('y', 'yes')
    return True if choice.strip().lower() in values else False

def floor(myvector, binsize):
    # rather than floor to the next lowest integer (i.e. multiple of 1), floor to the next lowest multiple of binsize
    return np.floor(myvector / binsize) * binsize

def ceil(myvector, binsize):
    # rather than ceil to the next highest integer (i.e. multiple of 1), ceil to the next highest multiple of binsize
    return np.ceil(myvector / binsize) * binsize

