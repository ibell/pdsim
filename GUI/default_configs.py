# -*- coding: utf-8 -*-
from __future__ import print_function

from families.scroll import scroll_yaml
from families.recip_compressor import recip_yaml
import yaml

def get_defaults(family):
    if family == 'scroll':
        return yaml.load(scroll_yaml,Loader=yaml.FullLoader)
    elif family == 'recip':
        return yaml.load(recip_yaml,Loader=yaml.FullLoader)
    else:
        raise ValueError('Your machine family [{f:s}] was not found'.format(f=family))

if __name__=='__main__':
    print(get_defaults('recip'))
