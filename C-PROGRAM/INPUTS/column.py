#!/usr/bin/env python3

import yaml

P = {'P': 0.00, 'Units': 'MPa', 'Definition': 'operating pressure'}

Dc = {'Dc': 0.00, 'Units': 'm', 'Definition': 'column diameter'}

drec = {'drec': 0.00, 'Units': 'm', 'Definition': 'recycle line diameter'}

vsep = {'vsep': 0.00, 'Units': 'm^3', 'Definition': 'gas-liquid separator volume'}

hset = {'hset': 0.00, 'Units': 'm', 'Definition': 'desired bed height'}

parameters = [P, Dc, drec, vsep, hset]

with open('column.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
