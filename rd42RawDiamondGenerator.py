#!/usr/bin/env python
import os, sys, shutil
sys.path.append('/home/sandiego/rd42AnalysisTools')  # TODO: HARDCODED!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from ConfigParser import ConfigParser
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import scipy.stats as sps
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
# from NoiseExtraction import NoiseExtraction
from Utils import *
import subprocess as subp
import multiprocessing as mp
import time

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12}
ped_axis = {'min': 0, 'max': 2**12}
cm_axis = {'min': -100, 'max': 100}


class RD42RawDiamondGenerator:
	def __init__(self, testnum, chs=10, events=1000000, landau_pos=100., landau_sc=20., noise=8., cm=9., buff=500):
		print 'Creating RD42 Raw Data'
		self.test_num = testnum
		self.ped_mean = 1000
		self.ped_spread = 100
		self.event_probab = 0.1
		self.chs = chs
		self.ped_buffer = buff
		self.maxEvents = events + self.ped_buffer
		self.landau_pos = landau_pos
		self.landau_sc = landau_sc
		self.noise = noise
		self.cm = cm
		self.raw_file = ro.TFile()
		self.raw_tree = ro.TTree()
		self.evt = np.zeros(1, 'uint32')
		self.adcs = np.zeros(self.chs, 'uint16')
		self.peds = np.zeros(self.chs, 'uint16')
		self.bar = None
		self.time = time.time()

	def CreateRawFile(self):
		self.raw_file = ro.TFile('testData{n}.root'.format(n=self.test_num), 'recreate')
		self.raw_tree = ro.TTree('rawTree', 'rawTree')
		self.raw_tree.Branch('event', self.evt, 'event/i')
		self.raw_tree.Branch('adc', self.adcs, 'adc[{n}]/s'.format(n=self.chs))

		self.peds = np.floor(0.5 + sps.norm.rvs(self.ped_mean, self.ped_spread, self.chs)).astype('uint16')
		print 'Filling tree "rawTree" in file testData{n}.root:'.format(n=self.test_num)
		self.bar = CreateProgressBarUtils(self.maxEvents)
		self.bar.start()
		for ev in xrange(self.maxEvents):
			self.evt.fill(ev)
			cmev = sps.norm.rvs(0, self.cm)
			noiseev = sps.norm.rvs(0, self.noise, self.chs)
			adcev = [self.peds[it] + noiseev[it] + cmev for it in xrange(self.chs)]
			if ev >= self.ped_buffer:
				hit = (sps.uniform.rvs(0, 1) <= self.event_probab)
				hitch = int(np.floor(sps.uniform.rvs(0, self.chs)))
				if hit:
					adcev[hitch] += sps.moyal.rvs(self.landau_pos, self.landau_sc)
					adcev[hitch] = 4095 if adcev[hitch] > 4095 else adcev[hitch]
			np.putmask(self.adcs, np.ones(self.chs, '?'), np.array(adcev, 'uint16'))
			self.raw_tree.Fill()
			self.bar.update(ev + 1)
		self.bar.finish()

		self.raw_file.Write()
		self.raw_file.Close()

		print 'Finished filling the rawTree in', time.time() - self.time, 'seconds :)'


def main():
	parser = OptionParser()
	parser.add_option('-n', '--testnumber', dest='testnum', type='int', help='test number')
	parser.add_option('-c', '--channels', dest='chs', type='int', default=10, help='Number of channels on the detector')
	parser.add_option('-e', '--events', dest='evts', type='int', default=1000000, help='Number of events to simulate')
	parser.add_option('-p', '--landaupos', dest='landaupos', type='float', default=100, help='position of the landauish distribution (moyal dist location parameter)')
	parser.add_option('-s', '--landauscale', dest='landauscale', type='float', default=20, help='landauish distribution width (moyal dist scale parameter)')
	parser.add_option('--noise', dest='noise', type='float', default=8, help='inherent noise of each channel in adcs')
	parser.add_option('--cm', dest='cm', type='float', default=9, help='common mode distribution that effects each event (in adcs)')
	parser.add_option('--buffer', dest='buffer', type='int', default=500, help='initial events without hits, just noise')

	(options, args) = parser.parse_args()
	testnum = int(options.testnum)
	chs = int(options.chs)
	events = int(options.evts)
	landau_pos = float(options.landaupos)
	landau_sc = float(options.landauscale)
	noise = float(options.noise)
	cm = float(options.cm)
	buff = int(options.buffer)

	rawGen = RD42RawDiamondGenerator(testnum=testnum, chs=chs, events=events, landau_pos=landau_pos, landau_sc=landau_sc, noise=noise, cm=cm, buff=buff)
	return rawGen

if __name__ == '__main__':
	rawGen = main()
	rawGen.CreateRawFile()
