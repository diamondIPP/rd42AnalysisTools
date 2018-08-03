#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
# from numpy import array, floor, average, std
import glob
import ROOT as ro
import numpy as np
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
		self.workingDir = os.getcwd()
		ro.gStyle.SetPalette(55)
		ro.gStyle.SetNumberContours(999)
		if not os.path.isdir('{d}/{r}/transparentAnalysis/root'.format(d=self.dir, r=self.run)):
			ExitMessage('The path {d}/{r}/transparentAnalysis/root does not exist. Exiting'.format(d=self.dir, r=self.run))
		os.chdir('{d}/{r}/transparentAnalysis/root'.format(d=self.dir, r=self.run))

	def ExtractHistoAndPlot(self, nameFile, nameDic, canvas_prefix='cRoot_', histoname='', option='colz'):
		areaFiles = glob.glob(nameFile)
		self.blaf[nameDic] = {}
		self.blac[nameDic] = {}
		self.blah[nameDic] = {}
		self.blac1[nameDic] = {}
		for areai, areaFile in enumerate(areaFiles):
			histoname = areaFile.split('.{r}'.format(r=self.run))[0] if histoname == '' else histoname
			self.blaf[nameDic][areai] = ro.TFile(areaFile)
			self.blac[nameDic][areai] = self.blaf[nameDic][areai].Get(canvas_prefix + areaFile.split('.{r}'.format(r=self.run))[0])
			self.blah[nameDic][areai] = self.blac[nameDic][areai].GetPrimitive(histoname)
			self.blac1[nameDic][areai] = ro.TCanvas(nameDic + '_area_' + str(areai), nameDic + '_area_' + str(areai), 1)
			self.blac1[nameDic][areai].cd()
			self.blah[nameDic][areai].Draw(option)

	def GetHighestChannelVsChargeSize1and2(self):
		clusterFiles = glob.glob('cNonHitPulseHeightDitribution_*OutOf10*')
		canvasNames = [clustf.split('.{r}'.format(r=self.run))[0] for clustf in clusterFiles]
		histoNames = [{'noCMC': 'h' + caNa[1:].replace('_',''), 'CMC': 'h' + caNa[1:].replace('_','').replace('tion','tionCMC')} for caNa in canvasNames]
		self.blaf['noise'] = {}
		self.blac['noise'] = {}
		self.blah['noise'] = {}
		self.blac1['noise'] = {}
		for i, filei in enumerate(clusterFiles):
			self.blaf['noise'][i] = ro.TFile(filei)
			self.blac['noise'][i] = self.blaf['noise'][i].Get(canvasNames[i])
			self.blah['noise'][i] = {}
			self.blac1['noise'][i] = {}
			for keyi, histoi in histoNames[i].iteritems():
				self.blah['noise'][i][keyi] = self.blac['noise'][i].GetPrimitive(histoi)
				self.blac1['noise'][i][keyi] = ro.TCanvas('noise_' + str(i) + '_' + keyi, 'noise_' + str(i) + '_' + keyi, 1)
				self.blah['noise'][i][keyi].Draw()



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
	bla.GetHighestChannelVsChargeSize1and2()

