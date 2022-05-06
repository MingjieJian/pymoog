#!/usr/bin/python
import os
import subprocess
import warnings
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
from shutil import copyfile
from tqdm.notebook import tqdm
from itertools import compress
from datetime import datetime

from scipy.spatial import Delaunay
from astropy import units as u
from astropy.modeling.models import BlackBody

def D2E(string):
    return string.replace('D', 'E')