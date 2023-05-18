#!/usr/bin/env python3

import yaml

absTol = {'absTol': 1e-5, 'Units': 'problem dependent', 'Definition': 'absolute tolerance'}

relTol = {'relTol': 1e-8, 'Units': 'non-dimensional', 'Definition': 'relative tolerance'}

maxIter = {'maxIter': 250, 'Units': 'non-dimensional', 'Definition': 'maximum number of iterations'}

parameters = [absTol, relTol, maxIter]

with open('convergence.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
