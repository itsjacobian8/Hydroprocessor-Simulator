#!/usr/bin/env python3

import yaml

sigma = {'sigma': 0.00, 'Units': 'N/m', 'Definition': 'gas-liquid interfacial tension'}

parameters = [sigma]

with open('interfacial.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
