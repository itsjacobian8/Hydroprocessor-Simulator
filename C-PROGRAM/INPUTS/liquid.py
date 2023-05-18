#!/usr/bin/env python3

import yaml

qlvr = {'qlvr': 0.0, 'Units': 'm^3/s', 'Definition': 'inlet liquid flow rate'}

rhol = {'rhol': 0.0, 'Units': 'kg/m^3', 'Definition': 'liquid density'}

mul = {'mul': 0.00, 'Units': 'Pa.s', 'Definition': 'liquid viscosity'}

parameters = [qlvr, rhol, mul]

with open('liquid.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
