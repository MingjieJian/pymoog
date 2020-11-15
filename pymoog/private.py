#!/usr/bin/python
import os
import subprocess
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import Delaunay
from astropy.modeling.blackbody import blackbody_lambda