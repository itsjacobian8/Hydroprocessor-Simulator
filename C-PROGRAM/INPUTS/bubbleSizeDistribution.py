#!/usr/bin/env python3

import yaml

classWidth = {'classWidth': 0.001, 'Units': 'mm', 'Definition': 'bubble class width'}

minCutOff = {'minCutOff': 0.025, 'Units': 'non-dimensional', 'Definition': 'minimum cumulative probability (for constructing bubble classes)'}

maxCutOff = {'maxCutOff': 0.975, 'Units': 'non-dimensional', 'Definition': 'maximum cumulative probability (for constructing bubble classes)'}

parameters = [classWidth, minCutOff, maxCutOff]

with open('bubbleSizeDistribution.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
