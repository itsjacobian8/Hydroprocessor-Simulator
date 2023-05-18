#!/usr/bin/env python3

import yaml

nr = {'nr': 0.00, 'Units': 'non-dimensional', 'Definition': 'number of risers'}

nor = {'nor': 0.0, 'Units': 'non-dimensional', 'Definition': 'number of outlet orifices per riser'}

dor = {'dor': 0.0, 'Units': 'm', 'Definition': 'distributor (i.e., bubble cap grid) outlet orifice diameter'}

parameters = [nr, nor, dor]

with open('grid.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
