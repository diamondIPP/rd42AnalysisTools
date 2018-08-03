#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
from collections import OrderedDict
from NoiseExtraction import NoiseExtraction
import os, sys, shutil
from Utils import *
import cPickle as pickle

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 100}
adc_axis = {'min': 0, 'max': 2**12}
ped_axis = {'min': 0, 'max': 2**12}
cm_axis = {'min': -100, 'max': 100}
ev_axis = {'min': 0, 'max': 0, 'bins': 1000}
ch_axis = {'min': -0.5, 'max': diaChs - 0.5, 'bins': diaChs}

ro.gStyle.SetPalette(55)
ro.gStyle.SetNumberContours(999)

class PHTransp:
	def __init__(self, run=22011, dir=''):
		self.dir = dir
		self.run = run
		self.blaf = {}
		self.blac = {}
		self.blah = {}
		self.blahp = {}
		self.blac1 = {}

	def GetXTransparentProfiles(self):
		self.GetTransparentProfiles('X')

	def GetYTransparentProfiles(self):
		self.GetTransparentProfiles('Y')

	def GetTransparentProfiles(self, var='X'):
		#ipdb.set_trace()
		self.blaf[var] = ro.TFile('{d}/{r}/transparentAnalysis/root/hLandau2HighestPredHit{v}_2outOf10.{r}.root'.format(d=self.dir, v=var, r=self.run))
		self.blac[var] = self.blaf[var].Get('cRoot_hLandau2HighestPredHit{v}_2outOf10'.format(v=var))
		self.blah[var] = self.blac[var].GetPrimitive('hLandau2HighestPredHit{v}_2outOf10'.format(v=var))
		self.blahp[var] = self.blah[var].ProfileY('PH2outOf10_{v}_mean'.format(v=var))
		self.blahp[var].SetTitle('PH2outOf10_{r}_{v}_mean'.format(r=self.run, v=var))
		self.blac1[var] = ro.TCanvas('blac1{v}'.format(v=var), 'blac1{v}'.format(v=var), 1)
		self.blac1[var].cd()
		self.blahp[var].Draw()

	def Get2DPHMap(self):
		pass


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', default=22022, type='int', help='Run to be analysed (e.g. 22022)')
	parser.add_option('-d', '--dir', dest='dir', default='.', type='string', help='source folder containing processed data of different runs')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

	bla = PHTransp(run, dir)
	bla.GetXTransparentProfiles()
	bla.GetYTransparentProfiles()
