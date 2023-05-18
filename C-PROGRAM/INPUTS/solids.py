#!/usr/bin/env python3

import yaml

mcat = {'mcat': 0.0, 'Units': 'kg', 'Definition': 'mass of catalyst inventory'}

rhop = {'rhop': 0.0, 'Units': 'kg/m^3', 'Definition': 'particle density'}

Vp = {'Vp': 0.00, 'Units': 'm^3', 'Definition': 'particle volume'}

Ap = {'Ap': 0.00, 'Units': 'm^2', 'Definition': 'particle surface area'}

parameters = [mcat, rhop, Vp, Ap]

with open('solids.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
