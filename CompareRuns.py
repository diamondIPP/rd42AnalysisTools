#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
import itertools as itt
from TransparentGrid import TransparentGrid
from TestAreas import TestAreas
from optparse import OptionParser
from Utils import *
from ConfigParser import ConfigParser

color_index = 0
hspace = 0.1
vspace = 0.0
class CompareRuns:
	def __init__(self, config_file=''):
		self.config_file = config_file
		# self.runlist = []
		self.outdir = '.'
		self.runs_path = {}
		self.runs_settings = {}
		self.run_numbers = []
		self.runs_subidr = {}
		# self.cellsize = cellsize
		# self.numstrips = numstrips
		self.testnum = 0
		self.ReadCompareConfig()
		self.subdir = 'test' + str(self.testnum)
		self.do_fit = False
		# self.ReadRunList()
		self.runs_ta = {}
		self.canvas = {}
		self.histo = {}
		self.stuff = []
		self.LoadRuns()
		self.run_colors = {run: ro.TColor(color_index + run, ReturnRGB(run, min(self.run_numbers), max(self.run_numbers))[0], ReturnRGB(run, min(self.run_numbers), max(self.run_numbers))[1], ReturnRGB(run, min(self.run_numbers), max(self.run_numbers))[2]) for run in self.run_numbers}
		self.numruns = len(self.run_numbers)
		self.run_pairs_comb = [[run1, run2] for run1, run2 in itt.combinations(self.run_numbers, 2)]
	#self.CompareAllPairs('ph1_{t}'.format(t=self.subdir))
	#self.CompareAllPairs('ph1_{t}'.format(t=self.subdir), do_norm=True)
	#self.CompareAllPairs('ph2_{t}'.format(t=self.subdir))
	#self.CompareAllPairs('ph2_{t}'.format(t=self.subdir), do_norm=True)
	#self.CompareAllPairs('PH_Saturated_Events')
	#self.CompareAllPairs('PH_Saturated_Events', do_norm=True)
	#self.SaveCanvasInList(self.canvas.keys())

	def ReadCompareConfig(self):
		if os.path.isfile(self.config_file):
			pars = ConfigParser()
			pars.read(self.config_file)
			print 'Reading config file for comparing runs...', ; sys.stdout.flush()
			if pars.has_section('RUNS'):
				if pars.has_option('RUNS', 'rundirs'):
					temp_dirs = pars.get('RUNS', 'rundirs')
					elements = temp_dirs.replace('{', '').replace('}', '')
					elements = elements.split(',')
					temp_runs_path = [Correct_Path(element) for element in elements if os.path.isdir(Correct_Path(element)) and (len(element) > 2)]
					self.run_numbers = np.array([int(runpath.split('/')[-1]) for runpath in temp_runs_path], 'uint32')
					self.runs_path = {self.run_numbers[i]: temp_runs_path[i] for i in xrange(len(temp_runs_path))}
					self.runs_subidr = {run: path.split('/' + str(run))[0] for run, path in self.runs_path.iteritems()}
				if pars.has_option('RUNS', 'settings'):
					temp_setts = pars.get('RUNS', 'settings')
					elements = temp_setts.replace('{', '').replace('}', '')
					elements = elements.split(',')
					self.runs_settings = [Correct_Path(element) for element in elements if os.path.isfile(Correct_Path(element)) and (len(element) > 2)]
				if pars.has_option('RUNS', 'outdir'):
					temp_outdir = pars.get('RUNS', 'outdir')
					self.outdir = Correct_Path(temp_outdir)
				if pars.has_option('RUNS', 'do_fit'):
					self.do_fit = pars.getboolean('RUNS', 'do_fit')
				if pars.has_option('RUNS', 'testnum'):
					self.testnum = pars.getint('RUNS', 'testnum')
			print 'Done'
	# def ReadRunList(self):
	# 	with open(self.runlist, 'r') as rl:
	# 		lines = rl.readlines()
	# 		self.runs_path = [Correct_Path(line.replace('\n', '')) for line in lines if os.path.isdir(Correct_Path(line.replace('\n', ''))) and ('#' not in line) and (';' not in line) and (len(line) > 2)]
	# 		self.run_numbers = np.array([int(runpath.split('/')[-1]) for runpath in self.runs_path], 'uint32')
	# 		self.runs_subidr = [path.split('/' + str(self.run_numbers[it]))[0] for it, path in enumerate(self.runs_path)]

	def LoadRuns(self):
		for run in self.run_numbers:
			if os.path.isdir(self.runs_path[run]):
				# self.runs_ta[run] = TestAreas(self.testnum, self.numstrips, self.runs_subidr[it], run, self.cellsize, self.do_fit)
				self.runs_ta[run] = TestAreas(self.runs_settings[0], run)
				self.runs_ta[run].trans_grid.SetLines()
				self.runs_ta[run].trans_grid.CreateTCutGs()
				self.runs_ta[run].SetTest()
				self.runs_ta[run].trans_grid.LoadPlotsInSubdir()

	def CompareAllPairs(self, histoname='ph2_test1', do_norm=False, plot_option='e hist'):
		for it, (run1, run2) in enumerate(self.run_pairs_comb):
			self.PlotPairs(run1, run2, histoname, do_norm, plot_option)

	def PlotPairs(self, run1, run2, histoname='ph2_test1', do_norm=False, plot_option='e hist'):
		tg1, tg2 = self.runs_ta[run1].trans_grid, self.runs_ta[run2].trans_grid
		suffix = '' if not do_norm else '_norm'
		if tg1.histo.has_key(histoname) and tg2.histo.has_key(histoname):
			new_name = histoname + '_{r1}_&_{r2}'.format(r1=run1, r2=run2) + suffix
			histo1name = histoname + '_' + str(run1)
			histo2name = histoname + '_' + str(run2)
			if not self.histo.has_key(histo1name):
				histo1 = tg1.histo[histoname]
				histo1.SetNameTitle(histo1name, histo1name)
				# pos1 = np.where(self.run_numbers == run1)[0][0]
				# histo1.FindObject('stats').SetX1NDC((hspace * (self.numruns - 2.0 * pos1) + pos1) / self.numruns)
				# histo1.FindObject('stats').SetX2NDC((hspace * (self.numruns - 2.0 * (pos1 + 1)) + pos1 + 1) / self.numruns)
				# histo1.FindObject('stats').SetY2NDC(1 - vspace)
				# histo1.FindObject('stats').SetY1NDC(1 - vspace - 0.3)
				self.histo[histo1name] = histo1
			if not self.histo.has_key(histo2name):
				histo2 = tg2.histo[histoname]
				histo2.SetNameTitle(histo2name, histo2name)
				# pos2 = np.where(self.run_numbers == run2)[0][0]
				# histo2.FindObject('stats').SetX1NDC((hspace * (self.numruns - 2.0 * pos2) + pos2) / self.numruns)
				# histo2.FindObject('stats').SetX2NDC((hspace * (self.numruns - 2.0 * (pos2 + 1)) + pos2 + 1) / self.numruns)
				# histo2.FindObject('stats').SetY2NDC(1 - vspace)
				# histo2.FindObject('stats').SetY1NDC(1 - vspace - 0.3)
				self.histo[histo2name] = histo2
			histo1 = self.histo[histo1name].Clone(histoname + '_{r1}_not{r2}'.format(r1=run1, r2=run2) + suffix)
			histo2 = self.histo[histo2name].Clone(histoname + '_{r2}_not{r1}'.format(r1=run1, r2=run2) + suffix)
			if do_norm:
				histo1.Sumw2(True)
				histo2.Sumw2(True)
				histo1.Scale(1.0 / histo1.Integral())
				histo2.Scale(1.0 / histo2.Integral())
				histo1.GetYaxis().SetTitle('norm. to 1')
				histo2.GetYaxis().SetTitle('norm. to 1')
			histo1.SetLineColor(self.run_colors[run1].GetNumber())
			histo1.SetMarkerColor(self.run_colors[run1].GetNumber())
			histo2.SetLineColor(self.run_colors[run2].GetNumber())
			histo2.SetMarkerColor(self.run_colors[run2].GetNumber())
			self.canvas[new_name] = ro.TCanvas('c_' + new_name, 'c_' + new_name, 1)
			self.canvas[new_name].SetGridx()
			self.canvas[new_name].SetGridy()
			self.canvas[new_name].SetTicky()
			ro.gStyle.SetOptTitle(0)
			hmax = histo1.GetBinContent(histo1.GetMaximumBin()) if histo1.GetBinContent(histo1.GetMaximumBin()) >= histo2.GetBinContent(histo2.GetMaximumBin()) else histo2.GetBinContent(histo2.GetMaximumBin())
			histo1.Draw(plot_option)
			histo2.Draw(plot_option + 'sames')
			# histo1.GetYaxis().SetRangeUser(0, 0.8 * hmax / ((1 - vspace - 0.3) - 0.1))
			histo1.GetYaxis().SetRangeUser(0, 0.8 * hmax / ((0.9) * 0.8))
			histo2.GetYaxis().SetRangeUser(0, 0.8 * hmax / ((0.9) * 0.8))
			histo2.FindObject('stats').SetY2NDC(0.6)
			histo2.FindObject('stats').SetY1NDC(0.3)
			ro.gPad.Update()
			legend = self.canvas[new_name].BuildLegend()
			ro.gPad.Update()
			legend.SetX1NDC(0.6)
			legend.SetX2NDC(0.9)
			legend.SetY1NDC(0.2)
			legend.SetY2NDC(0.3)
			ro.gPad.Update()
			self.stuff.append(histo1)
			self.stuff.append(histo2)

	def PlotSaturationEvents(self, histoname=''):
		tg_dic = {run: self.runs_ta[run].trans_grid for run in self.run_numbers}
		saturated = {}
		saturated_err = {}
		not_saturated = {}
		not_saturated_err = {}
		saturated_norm = {}
		saturated_err_norm = {}
		not_saturated_norm = {}
		not_saturated_err_norm = {}

		saturated_cuts = {}
		not_saturated_cuts = {}
		for run, tg in tg_dic.iteritems():
			num_strips_tree = int(tg.trans_tree.GetMaximum('numStrips'))
			satChs = '(' + '||'.join(['(diaChADC[clusterChannel{i}]==4095)'.format(i=i) for i in xrange(num_strips_tree)]) + ')'
			notSatChs = '(' + '&&'.join(['(diaChADC[clusterChannel{i}]<4095)'.format(i=i) for i in xrange(num_strips_tree)]) + ')'
			areaCut = '({c})'.format(c=tg.gridAreas.goodAreasCutNames_simplified_diamond)
			listSat = ['(transparentEvent)', satChs, areaCut]
			listNoSat = ['(transparentEvent)', notSatChs, areaCut]

			saturated[run] = tg.trans_tree.Draw('clusterChargeN', '(' + '&&'.join(listSat) + ')', 'goff')
			saturated_err[run] = np.sqrt(saturated[run])
			not_saturated[run] = tg.trans_tree.Draw('clusterChargeN', '(' + '&&'.join(listNoSat) + ')', 'goff')
			not_saturated_err[run] = np.sqrt(not_saturated[run])

			saturated_norm[run] = np.divide(saturated[run], saturated[run] + not_saturated[run], dtype='f8')
			saturated_err_norm[run] = np.divide(saturated_err[run] * not_saturated[run], np.power(saturated[run] + not_saturated[run], 2, dtype='f8'), dtype='f8')
			not_saturated_norm[run] = np.divide(not_saturated[run], saturated[run] + not_saturated[run], dtype='f8')
			not_saturated_err_norm[run] = np.divide(not_saturated_err[run] * saturated[run], np.power(saturated[run] + not_saturated[run], 2, dtype='f8'), dtype='f8')
		xvals = np.array(self.run_numbers, 'f8')

		y_sat_values = np.array([saturated[run] for run in self.run_numbers], 'f8')
		y_sat_values_err = np.array([saturated_err[run] for run in self.run_numbers], 'f8')
		y_not_sat_values = np.array([not_saturated[run] for run in self.run_numbers], 'f8')
		y_not_sat_values_err = np.array([not_saturated_err[run] for run in self.run_numbers], 'f8')

		y_sat_values_norm = np.array([saturated_norm[run] for run in self.run_numbers], 'f8')
		y_sat_values_err_norm = np.array([saturated_err_norm[run] for run in self.run_numbers], 'f8')
		y_not_sat_values_norm = np.array([not_saturated_norm[run] for run in self.run_numbers], 'f8')
		y_not_sat_values_err_norm = np.array([not_saturated_err_norm[run] for run in self.run_numbers], 'f8')

		graph_sat = ro.TGraphErrors(len(self.run_numbers), xvals, y_sat_values, np.zeros(len(self.run_numbers), 'f8'), y_sat_values_err)
		satname = 'SaturatedAfterCuts'
		graph_sat.SetNameTitle(satname, satname)
		graph_sat.GetXaxis().SetTitle('run')
		graph_sat.GetYaxis().SetTitle('# events')
		graph_sat.SetMarkerStyle(7)
		graph_sat.SetMarkerColor(ro.kRed)
		graph_nosat = ro.TGraphErrors(len(self.run_numbers), xvals, y_not_sat_values, np.zeros(len(self.run_numbers), 'f8'), y_not_sat_values_err)
		notsatname = 'NotSaturatedAfterCuts'
		graph_nosat.SetNameTitle(notsatname, notsatname)
		graph_nosat.GetXaxis().SetTitle('run')
		graph_nosat.GetYaxis().SetTitle('# events')
		graph_nosat.SetMarkerStyle(7)
		graph_nosat.SetMarkerColor(ro.kBlue)
		graph_sat_norm = ro.TGraphErrors(len(self.run_numbers), xvals, y_sat_values_norm, np.zeros(len(self.run_numbers), 'f8'), y_sat_values_err_norm)
		satnamenorm = 'SaturatedAfterCutsNorm'
		graph_sat_norm.SetNameTitle(satnamenorm, satnamenorm)
		graph_sat_norm.GetXaxis().SetTitle('run')
		graph_sat_norm.GetYaxis().SetTitle('norm')
		graph_sat_norm.GetYaxis().SetRangeUser(0, 1)
		graph_sat_norm.SetMarkerStyle(7)
		graph_sat_norm.SetMarkerColor(ro.kRed)
		graph_nosat_norm = ro.TGraphErrors(len(self.run_numbers), xvals, y_not_sat_values_norm, np.zeros(len(self.run_numbers), 'f8'), y_not_sat_values_err_norm)
		notsatnamenorm = 'NotSaturatedAfterCutsNorm'
		graph_nosat_norm.SetNameTitle(notsatnamenorm, notsatnamenorm)
		graph_nosat_norm.GetXaxis().SetTitle('run')
		graph_nosat_norm.GetYaxis().SetTitle('norm')
		graph_nosat_norm.GetYaxis().SetRangeUser(0, 1)
		graph_nosat_norm.SetMarkerStyle(7)
		graph_nosat_norm.SetMarkerColor(ro.kBlue)

		self.canvas[satname] = ro.TCanvas('c_' + satname, 'c_' + satname, 1)
		self.canvas[satname].SetGridx()
		self.canvas[satname].SetGridy()
		self.canvas[satname].SetTicky()
		graph_sat.Draw('ALP')

		self.canvas[notsatname] = ro.TCanvas('c_' + notsatname, 'c_' + notsatname, 1)
		self.canvas[notsatname].SetGridx()
		self.canvas[notsatname].SetGridy()
		self.canvas[notsatname].SetTicky()
		graph_nosat.Draw('ALP')

		self.canvas[satnamenorm] = ro.TCanvas('c_' + satnamenorm, 'c_' + satnamenorm, 1)
		self.canvas[satnamenorm].SetGridx()
		self.canvas[satnamenorm].SetGridy()
		self.canvas[satnamenorm].SetTicky()
		graph_sat_norm.Draw('ALP')

		self.canvas[notsatnamenorm] = ro.TCanvas('c_' + notsatnamenorm, 'c_' + notsatnamenorm, 1)
		self.canvas[notsatnamenorm].SetGridx()
		self.canvas[notsatnamenorm].SetGridy()
		self.canvas[notsatnamenorm].SetTicky()
		graph_nosat_norm.Draw('ALP')

		self.stuff.append(graph_sat)
		self.stuff.append(graph_nosat)
		self.stuff.append(graph_sat_norm)
		self.stuff.append(graph_nosat_norm)

	def SaveCanvasInList(self, list):
		if not os.path.isdir('{d}/{sd}'.format(d=self.outdir, sd=self.subdir)):
			os.makedirs('{d}/{sd}'.format(d=self.outdir, sd=self.subdir))
		for canvas in list:
			if self.canvas.has_key(canvas):
				self.canvas[canvas].SaveAs('{d}/{sd}/{c}.png'.format(d=self.outdir, sd=self.subdir, c=canvas))
				self.canvas[canvas].SaveAs('{d}/{sd}/{c}.root'.format(d=self.outdir, sd=self.subdir, c=canvas))

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-c', '--config', dest='config', type='string', help='Path to config file for comparing runs')
	# parser.add_option('-r', '--runlist', dest='runlist', type='string', help='Path to the file containing the list of runs')
	# parser.add_option('-o', '--outdir', dest='outdir', type='string', default='.' ,help='Path to the directory where the files will be saved')
	# parser.add_option('-c', '--cellsize', dest='cellsize', type='int', default=50, help='cell size of the square 3D device')
	# parser.add_option('-n', '--numstrips', dest='numstrips', type='int', default=2, help='Number of strips to use')
	# parser.add_option('-t', '--test', dest='testnumber', type='int', default=1, help='Run a automatically one of the predefined tests')
	# parser.add_option('-f', '--dofit', dest='dofit', default=False, action='store_true', help='Enables fitting')


	(options, args) = parser.parse_args()
	config = str(options.config)
	# runlist = str(options.runlist)
	# outdir = str(options.outdir)
	# cellsize = int(options.cellsize)
	# numstrips = int(options.numstrips)
	# testnum = int(options.testnumber)
	# do_fit = bool(options.dofit)

	c = CompareRuns(config)
