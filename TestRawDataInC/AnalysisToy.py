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
PedestalAnalysis = 'TestDataC'
inputInfoTypes = {'landau_pos': 'float', 'landau_sc': 'float', 'buff': 'int', 'ped_mean': 'int', 'ped_spread': 'float', 'noise': 'float', 'commonMode': 'float'}

class AnalysisToy:
	def __init__(self, filepath, threshold=5.):
		ro.gStyle.SetOptStat(1111)
		ro.gStyle.SetOptFit(111)
		print 'Analysing Toy raw data in', filepath
		self.filepath = Correct_Path(filepath)
		self.fileNameStem = self.filepath.split('/')[-1].split('.root')[0]
		self.fileNamePrefix = '/'.join(self.filepath.split('/')[:-1])
		self.threshold = threshold

		self.outFile = ro.TFile()
		self.boolIsOutFile = False
		self.outTree = ro.TTree()
		self.events = 0

		self.CheckOutFileExist()

		self.time = time.time()
		self.subp = None

		self.inputFileInfo = {'landau_pos': 0., 'landau_sc': 0., 'buff': 0, 'ped_mean': 0, 'ped_spread': 0., 'noise': 0., 'commonMode': 0.}

		self.canvas = {}
		self.histo = {}
		self.fits = {}
		self.w, self.window_shift = 0, 5

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self, canvas_name, self.w, self.window_shift)

	def CheckOutFileExist(self):
		if os.path.isfile('{d}/{f}Ped.root'.format(d=self.fileNamePrefix, f=self.fileNameStem)):
			self.boolIsOutFile = True

	def ExecutePedestalsCalculation(self):
		with open(os.devnull, 'w') as FNULL:
			command = os.getcwd() + '/' + PedestalAnalysis
			self.subp = subp.Popen([command, self.filepath, str(self.threshold)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
			while self.subp.poll() is None:
				time.sleep(2)
			if self.subp.poll() == 0:
				print 'Finished calculating the pedestals and signals :)'
			else:
				print 'Could have failed the calculation of pedestals and signals... Obtained return code:', self.subp.poll()
			CloseSubprocess(self.subp, True, False)

	def CreatePlots(self):
		print 'Creating plots...', ; sys.stdout.flush()
		self.time = time.time()
		self.OpenPedestalFile()

		nameh = 'signalPedestal'
		varx, vary = 'event', 'signal'
		xbins, fact = 1000000, 100.
		while xbins >= 10000:
			fact *= 10.
			xbins = int(RoundInt(self.events / fact, 'uint32'))
		xmin, xmax = 0, int(fact * xbins)
		ybins, ymin, ymax = 201, -100.5, 100.5
		cuts = self.CutOrs(['isPed'], [''])
		namex, namey, namez = 'event', 'signal[ADC]', 'entries'
		self.Draw2DHisto(nameh, varx, vary, xbins, xmin, xmax, ybins, ymin, ymax, cuts, namex, namey, namez)

		nameh = 'signalPedestalCMC'
		varx, vary = 'event', 'signalCMC'
		xbins, fact = 10000000, 100.
		while xbins >= 10000:
			fact *= 10.
			xbins = int(RoundInt(self.events / fact, 'uint32'))
		xmin, xmax = 0, int(fact * xbins)
		ybins, ymin, ymax = 201, -100.5, 100.5
		cuts = self.CutOrs(['isPedCMC'], [''])
		namex, namey, namez = 'event', 'signal_cmc[ADC]', 'entries'
		self.Draw2DHisto(nameh, varx, vary, xbins, xmin, xmax, ybins, ymin, ymax, cuts, namex, namey, namez)

		self.inputFileInfo['noise'] = self.histo[nameh].GetRMS(2)

		nameh = 'commonMode'
		varx, vary = 'event', 'cm'
		xbins, fact = 10000000, 100.
		while xbins >= 10000:
			fact *= 10.
			xbins = int(RoundInt(self.events / fact, 'uint32'))
		xmin, xmax = 0, int(fact * xbins)
		ybins, ymin, ymax = 201, -100.5, 100.5
		cuts = self.CutOrs(['isPedCMC'], [''])
		namex, namey, namez = 'event', 'common_mode[ADC]', 'entries'
		self.Draw2DHisto(nameh, varx, vary, xbins, xmin, xmax, ybins, ymin, ymax, cuts, namex, namey, namez)

		self.inputFileInfo['commonMode'] = self.histo[nameh].GetRMS(2)

		nameh = 'signalPH'
		var = 'signal'
		xbins, xmin, xmax = 500, 0, 5000
		cuts = self.CutAnds(['isPed', 'isSaturated'], ['==0', '==0'])
		namex, namey = 'signal_ph[ADC]', 'entries'
		self.Draw1DHisto(nameh, var, xbins, xmin, xmax, cuts, namex, namey)
		self.FitLandauish(nameh)
		ro.gPad.Update()

		nameh = 'signalPHCMC'
		var = 'signalCMC'
		xbins, xmin, xmax = 500, 0, 5000
		cuts = self.CutAnds(['isPedCMC', 'isSaturated'], ['==0', '==0'])
		namex, namey = 'signal_cmc_ph[ADC]', 'entries'
		self.Draw1DHisto(nameh, var, xbins, xmin, xmax, cuts, namex, namey)
		self.FitLandauish(nameh)
		ro.gPad.Update()

		self.inputFileInfo['landau_pos'], self.inputFileInfo['landau_sc'] = self.fits[nameh].Parameter(1), self.fits[nameh].Parameter(2)

		nameh = 'signalCMC'
		var = 'signalCMC'
		xbins, xmin, xmax = 5000, -100, 4900
		cuts = self.CutAnds(['isSaturated'], ['==0'])
		namex, namey = 'signal_cmc[ADC]', 'entries'
		self.Draw1DHisto(nameh, var, xbins, xmin, xmax, cuts, namex, namey)
		self.canvas[nameh].SetLogy()
		ro.gPad.Update()
		self.FitLandauishAndPed(nameh)
		ro.gPad.Update()

		print 'Done in', time.time() - self.time, 'seconds'

	def CutOrs(self, vars=[], vals=[]):
		return self.CutOp(vars, vals, '||')

	def CutAnds(self, vars=[], vals=[]):
		return self.CutOp(vars, vals, '&&')

	def CutOp(self, vars=[], vals=[], op='&&'):
		tvar = []
		tval = []
		for var, val in zip(vars, vals):
			if self.outTree.FindBranch(var):
				tvar.append(var)
				tval.append(val)
		cuts = ['(' + var + val + ')' for var, val in zip(tvar, tval)]
		return op.join(cuts)

	def Draw2DHisto(self, nameh, varx, vary, xbins, xmin, xmax, ybins, ymin, ymax, cuts, namex, namey, namez='entries'):
		self.canvas[nameh] = ro.TCanvas('c_' + nameh, 'c_' + nameh, 1)
		self.histo[nameh] = ro.TH2D('h_' + nameh, 'h_' + nameh, xbins, xmin, xmax, ybins, ymin, ymax)
		self.histo[nameh].GetXaxis().SetTitle(namex)
		self.histo[nameh].GetYaxis().SetTitle(namey)
		self.histo[nameh].GetZaxis().SetTitle(namez)
		self.outTree.Draw('{y}:{x}>>h_'.format(x=varx, y=vary) + nameh, cuts, 'colz')
		if self.histo[nameh].FindObject('stats'):
			self.histo[nameh].FindObject('stats').SeOptStat(2211)
		self.PosCanvas(nameh)

	def Draw1DHisto(self, nameh, var, xbins, xmin, xmax, cuts, namex, namey='entries'):
		self.canvas[nameh] = ro.TCanvas('c_' + nameh, 'c_' + nameh, 1)
		self.histo[nameh] = ro.TH1D('h_' + nameh, 'h_' + nameh, xbins, xmin, xmax)
		self.histo[nameh].GetXaxis().SetTitle(namex)
		self.histo[nameh].GetYaxis().SetTitle(namey)
		self.outTree.Draw('{v}>>h_'.format(v=var) + nameh, cuts, 'e hist')
		SetDefault1DCanvasSettings(self.canvas[nameh])
		ro.gPad.Update()
		SetDefault1DStats(self.histo[nameh])
		self.PosCanvas(nameh)
		ro.gPad.Update()

	def OpenPedestalFile(self):
		self.outFile = ro.TFile('{d}/{f}Ped.root'.format(d=self.fileNamePrefix, f=self.fileNameStem), 'read')
		self.outTree = self.outFile.Get('pedTree')
		self.events = self.outTree.GetEntries()

		for key in self.inputFileInfo.keys():
			if self.outFile.Get(key):
				self.inputFileInfo[key] = int(self.outFile.Get(key).GetTitle()) if inputInfoTypes[key] == 'int' else float(self.outFile.Get(key).GetTitle())

	def FitLandauish(self, name, iterations=2):
		"""
		Function for ROOT of a Landauish function based on moyal distribution
		"""
		if name in self.histo.keys():
			if self.histo[name]:
				if self.histo[name].GetEntries() > 0:
					if name in self.canvas.keys():
						self.canvas[name].cd()
					xmean, xrms, entries = self.histo[name].GetMean(), self.histo[name].GetRMS(), self.histo[name].GetEntries()
					if self.inputFileInfo['landau_pos'] == 0. or self.inputFileInfo['landau_sc'] == 0.:
						xmin, xmax = max(0, xmean - 3 * xrms), xmean + 5 * xrms
						loc, sc = max(1, xmean - xrms/2.), xrms / 2.
					else:
						loc, sc = self.inputFileInfo['landau_pos'], self.inputFileInfo['landau_sc']
						toyExp = sps.moyal.rvs(loc, sc, int(1e6))
						xmin, xmax = max(0, toyExp.min() - xrms), toyExp.max() + xrms
					func = ro.TF1('f_landauish_' + name, '[0]*TMath::Exp(-(x-[1])/(2*[2])-0.5*TMath::Exp(-(x-[1])/[2]))', 0, self.histo[name].GetBinLowEdge(1 + self.histo[name].FindLastBinAbove(0, 1)))
					func.SetNpx(1000)
					func.SetLineStyle(2)
					func.SetLineColor(ro.kRed)
					func.SetLineWidth(2)
					params = np.array([entries, loc, sc], 'float64')
					func.SetParameters(params)
					func.SetParNames('gain', 'pos', 'scale')
					for it in xrange(iterations):
						temp = self.histo[name].Fit('f_landauish_' + name, 'QIEBMSN', 'goff', xmin, xmax)
						params = np.array((temp.Parameter(0), temp.Parameter(1), temp.Parameter(2)), 'float64')
						func.SetParameters(params)
					self.fits[name] = self.histo[name].Fit('f_landauish_' + name, 'QIEBMS', 'goff', xmin, xmax)
					self.histo[name].GetFunction('f_landauish_' + name).Draw('same')
					ro.gPad.Update()

	def FitLandauishAndPed(self, name, iterations=2):
		"""
		Function for ROOT of a Landauish function based on moyal distribution
		"""
		if name in self.histo.keys():
			if self.histo[name]:
				if self.histo[name].GetEntries() > 0:
					if name in self.canvas.keys():
						self.canvas[name].cd()
					xmean, xrms, entries = self.histo[name].GetMean(), self.histo[name].GetRMS(), self.histo[name].GetEntries()
					(loc, sc) = (1000., 100.) if self.inputFileInfo['landau_pos'] == 0. or self.inputFileInfo['landau_sc'] == 0. else (self.inputFileInfo['landau_pos'], self.inputFileInfo['landau_sc'])
					sigm = 10. if self.inputFileInfo['noise'] == 0. else self.inputFileInfo['noise']
					xmin, xmax = self.histo[name].GetBinLowEdge(self.histo[name].FindFirstBinAbove(0, 1)), self.histo[name].GetBinLowEdge(1 + self.histo[name].FindLastBinAbove(0, 1))
					func = ro.TF1('f_ped_landauish_' + name, '[0]*TMath::Exp(-(x-[1])/(2*[2])-0.5*TMath::Exp(-(x-[1])/[2]))+gaus(3)', xmin, xmax)
					func.SetNpx(1000)
					func.SetLineStyle(2)
					func.SetLineColor(ro.kRed)
					func.SetLineWidth(2)
					params = np.array([entries/10., loc, sc, entries, 0., sigm], 'float64')
					func.SetParameters(params)
					func.SetParNames('Lgain', 'Lpos', 'Lscale', 'Ggain', 'Gpos', 'Gsigma')
					for it in xrange(iterations):
						temp = self.histo[name].Fit('f_ped_landauish_' + name, 'QIEBMSN', 'goff', xmin, xmax)
						params = np.array((temp.Parameter(0), temp.Parameter(1), temp.Parameter(2), temp.Parameter(3), temp.Parameter(4), temp.Parameter(5)), 'float64')
						func.SetParameters(params)
					self.fits[name] = self.histo[name].Fit('f_ped_landauish_' + name, 'QIEBMS', 'goff', xmin, xmax)
					self.histo[name].GetFunction('f_ped_landauish_' + name).Draw('same')
					ro.gPad.Update()

	def SavePlots(self):
		if not os.path.isdir('{d}/{p}'.format(d=self.fileNamePrefix, p=self.fileNameStem)):
			os.makedirs('{d}/{p}'.format(d=self.fileNamePrefix, p=self.fileNameStem))
		for canvas in self.canvas.iterkeys():
			if self.canvas.has_key(canvas):
				if self.canvas[canvas]:
					self.canvas[canvas].SaveAs('{d}/{p}/{c}.png'.format(d=self.fileNamePrefix, p=self.fileNameStem, c=canvas))
					self.canvas[canvas].SaveAs('{d}/{p}/{c}.root'.format(d=self.fileNamePrefix, p=self.fileNameStem, c=canvas))


def main():
	parser = OptionParser()
	parser.add_option('-i', '--input', dest='input', type='string', help='input file with raw data')
	parser.add_option('-t', '--threshold', dest='threshold', type='float', default=5, help='threshold in sigmas for pedestal determination')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true', help='toggles overwrite output file if it already exists')
	parser.add_option('-b', '--batch', dest='batch', default=False, action='store_true', help='enables batch mode (no print)')

	(options, args) = parser.parse_args()
	infile = str(options.input)
	threshold = float(options.threshold)
	force = bool(options.force)
	bat = bool(options.batch)

	if bat:
		ro.gROOT.SetBatch(ro.kTRUE)

	anaT = AnalysisToy(infile, threshold)
	if force or not anaT.boolIsOutFile:
		anaT.ExecutePedestalsCalculation()
	return anaT

if __name__ == '__main__':
	anaT = main()
	anaT.CreatePlots()
	anaT.SavePlots()
