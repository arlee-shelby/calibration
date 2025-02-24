import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from models.base import Session
from models.pixel import Pixel
import numpy as np

session = Session()

class PixelController:
    def get_pixels():
        return session.query(Pixel).all()

    def get_pixel_by_number(num):
        return session.query(Pixel).filter_by(pixel_number = num).all()[0]

    def energy_hist(pixel, bins):
        hist, bin_edges = np.histogram(pixel.energy, bins = bins)
        width = bin_edges[1]-bin_edges[0]
        return hist, bin_edges, width
        