#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#       
# Source class
# 

class Source:
    def __init__(self, name, asp_size, kc, Mw ):
        self.name = name
        self.asp_size = asp_size
        self.kc = kc
        self.Mw = Mw

    def __str__(self):
        source = " Source's name: " + str(self.name) + "\n"
        source += " asp_size (km): " + str(self.asp_size)  + "\n"
        source += " kc: " + str(self.asp_size) + " Mw: " + str(self.Mw) + "\n"




