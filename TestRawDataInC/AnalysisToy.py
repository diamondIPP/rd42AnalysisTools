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
		self.OpenPedestalFile()

		self.canvas['signalPedestal'] = ro.TCanvas('c_signalPedestal', 'c_signalPedestal', 1)
		self.histo['signalPedestal'] = ro.TH2D('h_signalPedestal', 'h_signalPedestal', int(RoundInt(self.events / 1000., 'uint32')), 0, int(1000 * RoundInt(self.events / 1000., 'uint32')), 201, -100.5, 100.5)
		self.histo['signalPedestal'].GetXaxis().SetTitle('event')
		self.histo['signalPedestal'].GetYaxis().SetTitle('signal[ADC]')
		self.histo['signalPedestal'].GetZaxis().SetTitle('entries')
		self.outTree.Draw('signal:event>>h_signalPedestal', 'isPed', 'colz')
		if self.histo['signalPedestal'].FindObject('stats'):
			self.histo['signalPedestal'].FindObject('stats').SeOptStat(2211)
		self.PosCanvas('signalPedestal')

		self.canvas['signalPedestalCMC'] = ro.TCanvas('c_signalPedestalCMC', 'c_signalPedestalCMC', 1)
		self.histo['signalPedestalCMC'] = ro.TH2D('h_signalPedestalCMC', 'h_signalPedestalCMC', int(RoundInt(self.events / 1000., 'uint32')), 0, int(1000 * RoundInt(self.events / 1000., 'uint32')), 201, -100.5, 100.5)
		self.histo['signalPedestalCMC'].GetXaxis().SetTitle('event')
		self.histo['signalPedestalCMC'].GetYaxis().SetTitle('signal CMC[ADC]')
		self.histo['signalPedestalCMC'].GetZaxis().SetTitle('entries')
		self.outTree.Draw('signalCMC:event>>h_signalPedestalCMC', 'isPedCMC', 'colz')
		if self.histo['signalPedestalCMC'].FindObject('stats'):
			self.histo['signalPedestalCMC'].FindObject('stats').SeOptStat(2211)
		self.PosCanvas('signalPedestalCMC')

		self.canvas['commonMode'] = ro.TCanvas('c_commonMode', 'c_commonMode', 1)
		self.histo['commonMode'] = ro.TH2D('h_commonMode', 'h_commonMode', int(RoundInt(self.events / 1000., 'uint32')), 0, int(1000 * RoundInt(self.events / 1000., 'uint32')), 201, -100.5, 100.5)
		self.histo['commonMode'].GetXaxis().SetTitle('event')
		self.histo['commonMode'].GetYaxis().SetTitle('common mode [ADC]')
		self.histo['commonMode'].GetZaxis().SetTitle('entries')
		self.outTree.Draw('cm:event>>h_commonMode', 'isPedCMC', 'colz')
		if self.histo['commonMode'].FindObject('stats'):
			self.histo['commonMode'].FindObject('stats').SeOptStat(2211)
		self.PosCanvas('commonMode')

		self.canvas['signalPH'] = ro.TCanvas('c_signalPH', 'c_signalPH', 1)
		self.histo['signalPH'] = ro.TH1D('h_signalPH', 'h_signalPH', 500, 0, 5000)
		self.histo['signalPH'].GetXaxis().SetTitle('signalPH[ADC]')
		self.histo['signalPH'].GetYaxis().SetTitle('entries')
		self.outTree.Draw('signal>>h_signalPH', '(isPed==0)&&(isSaturated==0)', 'e hist')
		# self.histo['signalPH'].FindObject('stats').SetOptFit(1)
		SetDefault1DCanvasSettings(self.canvas['signalPH'])
		ro.gPad.Update()
		SetDefault1DStats(self.histo['signalPH'])
		self.FitLandauish('signalPH')
		ro.gPad.Update()
		self.PosCanvas('signalPH')

		self.canvas['signalPHCMC'] = ro.TCanvas('c_signalPHCMC', 'c_signalPHCMC', 1)
		self.histo['signalPHCMC'] = ro.TH1D('h_signalPHCMC', 'h_signalPHCMC', 500, 0, 5000)
		self.histo['signalPHCMC'].GetXaxis().SetTitle('signalPHCMC[ADC]')
		self.histo['signalPHCMC'].GetYaxis().SetTitle('entries')
		self.outTree.Draw('signalCMC>>h_signalPHCMC', '(isPedCMC==0)&&(isSaturated==0)', 'e hist')
		# self.histo['signalPHCMC'].FindObject('stats').SetOptFit(1)
		SetDefault1DCanvasSettings(self.canvas['signalPHCMC'])
		ro.gPad.Update()
		SetDefault1DStats(self.histo['signalPHCMC'])
		self.FitLandauish('signalPHCMC')
		ro.gPad.Update()
		self.PosCanvas('signalPHCMC')

		self.canvas['signalCMC'] = ro.TCanvas('c_signalCMC', 'c_signalCMC', 1)
		self.histo['signalCMC'] = ro.TH1D('h_signalCMC', 'h_signalCMC', 5000, -100, 4900)
		self.histo['signalCMC'].GetXaxis().SetTitle('signalCMC[ADC]')
		self.histo['signalCMC'].GetYaxis().SetTitle('entries')
		self.outTree.Draw('signalCMC>>h_signalCMC', 'isSaturated==0', 'e hist')
		# self.histo['signalCMC'].FindObject('stats').SetOptFit(1)
		SetDefault1DCanvasSettings(self.canvas['signalCMC'])
		self.canvas['signalCMC'].SetLogy()
		ro.gPad.Update()
		SetDefault1DStats(self.histo['signalCMC'])
		self.FitLandauishAndPed('signalCMC')
		ro.gPad.Update()
		self.PosCanvas('signalCMC')
		print 'Done'

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
			if name in self.canvas.keys():
				self.canvas[name].cd()
			xmean, xrms, entries = self.histo[name].GetMean(), self.histo[name].GetRMS(), self.histo[name].GetEntries()
			if self.inputFileInfo['landau_pos'] == 0. or self.inputFileInfo['landau_sc'] == 0.:
				xmin, xmax = max(0, xmean - 3 * xrms), xmean + 3 * xrms
				loc, sc = max(1, xmean - xrms/2.), xrms / 2.
			else:
				loc, sc = self.inputFileInfo['landau_pos'], self.inputFileInfo['landau_sc']
				toyExp = sps.moyal.rvs(loc, sc, int(1e6))
				xmin, xmax = max(0, toyExp.min() - xrms), toyExp.max() + xrms
			func = ro.TF1('f_landauish_' + name, '[0]*TMath::Exp(-(x-[1])/(2*[2])-0.5*TMath::Exp(-(x-[1])/[2]))', 0, xmax)
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
			if name in self.canvas.keys():
				self.canvas[name].cd()
			xmean, xrms, entries = self.histo[name].GetMean(), self.histo[name].GetRMS(), self.histo[name].GetEntries()
			if self.inputFileInfo['landau_pos'] == 0. or self.inputFileInfo['landau_sc'] == 0.:
				loc, sc = 1000., 100.
				sigm = 10.
			else:
				loc, sc = self.inputFileInfo['landau_pos'], self.inputFileInfo['landau_sc']
				sigm = self.inputFileInfo['noise']
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
					# self.canvas[canvas].Print('{d}/{r}/{sd}/{c}.png'.format(d=self.dir, r=self.run, sd=self.pkl_sbdir, c=canvas))
					self.canvas[canvas].SaveAs('{d}/{p}/{c}.root'.format(d=self.fileNamePrefix, p=self.fileNameStem, c=canvas))


def main():
	parser = OptionParser()
	parser.add_option('-i', '--input', dest='input', type='string', help='input file with raw data')
	parser.add_option('-t', '--threshold', dest='threshold', type='float', default=5, help='threshold in sigmas for pedestal determination')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true', help='toggles overwrite output file if it already exists')

	(options, args) = parser.parse_args()
	infile = str(options.input)
	threshold = float(options.threshold)
	force = bool(options.force)

	anaT = AnalysisToy(infile, threshold)
	if force or not anaT.boolIsOutFile:
		anaT.ExecutePedestalsCalculation()
	return anaT

if __name__ == '__main__':
	anaT = main()
	anaT.CreatePlots()
	anaT.SavePlots()
