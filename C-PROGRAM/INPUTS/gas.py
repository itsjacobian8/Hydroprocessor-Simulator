#!/usr/bin/env python3

import yaml

qgtMin = {'qgtMin':0.00, 'Units':'m^3/s', 'Definition': 'minimum inlet gas flow rate'}

qgtMax = {'qgtMax':0.00, 'Units':'m^3/s', 'Definition': 'maximum inlet gas flow rate'}

qgtIncrement = {'qgtIncrement':0.00, 'Units':'m^3/s', 'Definition': 'step size for inlet gas flow rate'}

rhog = {'rhog': 0.00, 'Units': 'kg/m^3', 'Definition': 'gas density'}

parameters = [qgtMin, qgtMax, qgtIncrement, rhog]

with open('gas.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
