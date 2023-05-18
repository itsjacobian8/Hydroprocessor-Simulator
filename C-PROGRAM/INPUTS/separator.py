#!/usr/bin/env python3

import yaml


beta1 = {'beta1': 0.00, 'Units': '1/s', 'Definition': 'geometry dependent parameter of gas-liquid separation submodel'}

beta2 = {'beta2': 0.00, 'Units': '1/mm', 'Definition': 'bubble and geometry dependent parameter of gas-liquid separation submodel'}

Rstep = {'Rstep': 0.00, 'Units': 'non-dimensional', 'Definition': 'step size for the liquid recycle ratio'}

parameters = [beta1, beta2, Rstep]

with open('separator.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
