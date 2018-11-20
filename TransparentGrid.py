#!/usr/bin/env python
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
from copy import deepcopy
from collections import OrderedDict
import os, sys, shutil
from Utils import *
import cPickle as pickle

class TransparentGrid:
    def __init__(self, dir='', run=25209):
        ro.gStyle.SetPalette(55)
        ro.gStyle.SetNumberContours(999)
        ro.TFormula.SetMaxima(100000)
        self.run = run
        self.dir = os.path.abspath(os.path.expanduser(os.path.expandvars(dir)))
        self.trans_file = None
        self.trans_tree = None
        self.align_file = None
        self.align_obj = None
        self.align_info = {'xoff': float(0), 'phi': float(0)}
        self.num_cols = 19
        self.ch_ini = 0
        self.ch_end = 84
        self.phbins = 200
        self.phmin = 0
        self.phmax = 4000
        self.pkl = None
        self.loaded_pickle = False
        self.row_info = {'0': float(61.84669791829917), 'm': float(0.02551248435536136), 'num': 27, 'pitch': 1.008}
        self.vertical_lines = []
        self.vertical_lines_tline = []
        self.horizontal_lines = []
        self.horizontal_lines_tline = []
        self.bins_per_ch_x = 3
        self.bins_per_ch_y = 3
        self.canvas = {}
        self.profile = {}
        self.histo = {}
        self.names = []
        self.tcutgs = {}
        self.goodAreas = []
        self.badAreas = []
        self.goodAreasCutNames = ''
        self.badAreasCutNames = ''

    def CheckFoldersAndFiles(self):
        if not os.path.isdir(self.dir):
            ExitMessage('The path to the directory "{d}" does not exist. Exiting...'.format(d=self.dir), code=os.EX_NOINPUT)
        if not os.path.isdir('{d}/{r}'.format(d=self.dir, r=self.run)):
            ExitMessage('There is no run {r} in the directory "{d}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)
        if not os.path.isfile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run)):
            ExitMessage('There is no transparent root file "transparent.{r}.root" in the directory "{d}/{r}". Exiting...'.format(r=self.run, d=self.dir), code=os.EX_NOINPUT)

    def OpenFileAndGetTree(self):
        self.trans_file = ro.TFile('{d}/{r}/transparent.{r}.root'.format(d=self.dir, r=self.run), 'r')
        self.trans_tree = self.trans_file.Get('transparentTree')

    def FindDiamondChannelLimits(self):
        temph = ro.TH1F('temph', 'temph', 128, -0.5,127.5)
        self.trans_tree.Draw('diaChXPred>>temph', '', 'goff')
        self.ch_ini = int(temph.GetBinCenter(temph.FindFirstBinAbove(0, 1)))
        self.ch_end = int(temph.GetBinCenter(temph.FindLastBinAbove(0, 1)))
        self.num_cols = self.ch_end - self.ch_ini + 1

    def TryLoadPickle(self):
        picklepath = '{d}/{r}/y_low_params.{r}.pkl'.format(d=self.dir, r=self.run)
        if os.path.isfile(picklepath):
            with open(picklepath, 'rb') as pkl:
                self.pkl = pickle.load(pkl)
                self.loaded_pickle = True

    def SetLines(self, try_align=False):
        self.TryLoadPickle()
        if self.loaded_pickle:
            self.row_info = self.pkl['row_info']
            self.vertical_lines = self.pkl['vertical_lines']
            self.horizontal_lines = self.pkl['horizontal_lines']
            self.align_info = self.pkl['align_info']
        elif try_align:
            self.FindHorizontalParametersThroughAlignment()
            # self.CreateLines()
        else:
            self.AskUserLowerYLine()
            self.CreateLines()

    def AskUserLowerYLine(self):
        self.row_info['0'] = self.GetFromUser('Enter the y axis intercept in silicon space for the lower detector limit (scalar between 0 and 255): ', typ='float', limits=[0, 255])
        self.row_info['m'] = self.GetFromUser('Enter the slope for the lower detector line (scalar between -1 and 1): ', typ='float', limits=[-1, 1])
        self.row_info['num'] = self.GetFromUser('Enter the number of rows: ', typ='int', limits=[1, 1000])
        self.row_info['pitch'] = self.GetFromUser('Enter the effective pitch in sil space in Y axis: ', typ='float', limits=[0, 255])

    def GetFromUser(self, message, typ, limits=[]):
        cont = False
        tempv = 0
        while not cont:
            tempv = raw_input(message)
            if typ == 'int':
                if IsInt(tempv):
                    tempv = int(tempv)
                    if len(limits) == 2:
                        if limits[0] <= tempv <= limits[1]:
                            cont = True
                    else:
                        cont = True
            if typ == 'float':
                if IsFloat(tempv):
                    tempv = float(tempv)
                    if len(limits) == 2:
                        if limits[0] <= tempv <= limits[1]:
                            cont = True
                    else:
                        cont = True
        return tempv

    def FindHorizontalParametersThroughAlignment(self):
        self.LoadAlignmentParameters()

    def LoadAlignmentParameters(self):
        if os.path.isfile('{d}/{r}/alignment.{r}.root'.format(d=self.dir, r=self.run)):
            self.align_file = ro.TFile('{d}/{r}/alignment.{r}.root'.format(d=self.dir, r=self.run), 'r')
            self.align_obj = self.align_file.Get('alignment')
            self.align_info['xoff'] = self.align_obj.GetXOffset(4)
            self.align_info['phi'] = self.align_obj.GetPhiXOffset(4)

    def CreateLines(self):
        linev = self.GetVerticalLine(x=self.ch_ini - 0.5)
        lineh = self.GetHorizontalLine(y=self.row_info['0'] + self.row_info['m'] * (self.ch_ini - 0.5))
        self.vertical_lines.append(linev)
        self.horizontal_lines.append(lineh)
        for col in xrange(self.num_cols):
            linev = self.GetVerticalLine(self.ch_ini + col + 0.5)
            self.vertical_lines.append(linev)
        for row in xrange(self.row_info['num']):
            lineh = self.GetHorizontalLine(y=self.row_info['0'] + self.row_info['m'] * (self.ch_ini - 0.5) + (row + 1) * self.row_info['pitch'])
            self.horizontal_lines.append(lineh)

    def GetVerticalLine(self, x):
        return {0: {'x': x, 'y': self.row_info['0'] + self.row_info['m'] * x}, 1: {'x': x, 'y': self.row_info['0'] + self.row_info['m'] * x + self.row_info['num'] * self.row_info['pitch']}}

    def GetHorizontalLine(self, y):
        return {0: {'x': self.ch_ini - 0.5, 'y': y}, 1: {'x': self.ch_end + 0.5, 'y': y + self.row_info['m'] * (self.ch_end - self.ch_ini + 1)}}

    def CreateLinesTGraph(self):
        for lineh in self.horizontal_lines:
            self.horizontal_lines_tline.append(ro.TLine(lineh[0]['x'], lineh[0]['y'], lineh[1]['x'], lineh[1]['y']))
            self.horizontal_lines_tline[-1].SetLineColor(ro.kRed)
        for linev in self.vertical_lines:
            self.vertical_lines_tline.append(ro.TLine(linev[0]['x'], linev[0]['y'], linev[1]['x'], linev[1]['y']))
            self.vertical_lines_tline[-1].SetLineColor(ro.kRed)

    def DrawProfile2D(self, name, varz):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(128 * self.bins_per_ch_x), - 0.5, 127.5, int(256 * self.bins_per_ch_y), -0.5, 255.5)
        self.canvas[name] = ro.TCanvas('c_'+name, 'c_'+name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}:fidY:diaChXPred>>h_{n}'.format(z=varz, n=name), 'transparentEvent', 'colz prof')

    def DrawLines(self, name):
        self.canvas[name].cd()
        for lineh in self.horizontal_lines_tline:
            lineh.Draw('same')
        for linev in self.vertical_lines_tline:
            linev.Draw('same')

    def ResetLines(self):
        self.horizontal_lines = []
        self.horizontal_lines_tline = []
        self.vertical_lines = []
        self.vertical_lines_tline = []

    def CreateTCutGs(self):
        def GetNumpyArraysX(coli):
            x0 = self.ch_ini - 0.5 + coli
            x1 = self.ch_ini + 0.5 + coli
            return np.array((x0, x0, x1, x1, x0), 'f8')
        def GetNumpyArraysY(coli, rowi):
            y0 = self.row_info['0'] + self.row_info['m'] * (self.ch_ini - 0.5 + coli) + rowi * self.row_info['pitch']
            y1 = self.row_info['0'] + self.row_info['m'] * (self.ch_ini - 0.5 + coli) + (rowi + 1) * self.row_info['pitch']
            y2 = self.row_info['0'] + self.row_info['m'] * (self.ch_ini + 0.5 + coli) + (rowi + 1) * self.row_info['pitch']
            y3 = self.row_info['0'] + self.row_info['m'] * (self.ch_ini + 0.5 + coli) + rowi * self.row_info['pitch']
            return np.array((y0, y1, y2, y3, y0), 'f8')
        for col in xrange(self.num_cols):
            self.tcutgs[col] = {}
            for row in xrange(self.row_info['num']):
                tempx = GetNumpyArraysX(col)
                tempy = GetNumpyArraysY(col, row)
                self.tcutgs[col][row] = ro.TCutG('cutg_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
                self.tcutgs[col][row].SetNameTitle('cutg_{c}_{r}'.format(c=col, r=row), 'cutg_{c}_{r}'.format(c=col, r=row))
                self.tcutgs[col][row].SetVarX('diaChXPred')
                self.tcutgs[col][row].SetVarY('fidY')
                self.tcutgs[col][row].SetLineColor(ro.kRed)

    def AddGoodAreas(self, col, row):
        self.tcutgs[col][row].SetLineColor(ro.kRed)
        self.goodAreas.append(self.tcutgs[col][row])

    def AddBadAreas(self, col, row):
        self.tcutgs[col][row].SetLineColor(ro.kBlue)
        self.badAreas.append(self.tcutgs[col][row])

    def DrawGoodAreas(self, name):
        self.canvas[name].cd()
        for area in self.goodAreas:
            area.Draw('same')

    def DrawBadAreas(self, name):
        self.canvas[name].cd()
        for area in self.badAreas:
            area.Draw('same')

    def SelectGoodAndBadByThreshold(self, val=500):
        for col in xrange(self.num_cols):
            for row in xrange(self.row_info['num']):
                temph = ro.TH1F('temphrc', 'temphrc', 200, 0, 4000)
                self.trans_tree.Draw('clusterChargeN>>temphrc', 'transparentEvent&&({n})'.format(n=self.tcutgs[col][row].GetName()), 'goff')
                if temph.GetMean() > val:
                    self.AddGoodAreas(col, row)
                else:
                    self.AddBadAreas(col, row)
                del temph

        tempgood = [cut.GetName() for cut in self.goodAreas]
        self.goodAreasCutNames = '((' + ')||('.join(tempgood) + '))'
        tempbad = [cut.GetName() for cut in self.badAreas]
        self.badAreasCutNames = '((' + ')||('.join(tempbad) + '))'

    def ResetAreas(self):
        self.goodAreas = []
        self.badAreas = []

    def DrawPHGoodAreas(self, name, var='clusterChargeN'):
        self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, self.phbins, self.phmin, self.phmax)
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}>>h_{n}'.format(z=var, n=name), 'transparentEvent&&{n}'.format(n=self.goodAreasCutNames))

    def DrawPHBadAreas(self, name, var='clusterChargeN'):
        self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, self.phbins, self.phmin, self.phmax)
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}>>h_{n}'.format(z=var, n=name), 'transparentEvent&&{n}'.format(n=self.badAreasCutNames))

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
    parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
    parser.add_option('-a', '--al', dest='al', action='store_true', default=False, help='enable find grid through data and alignment')

    (options, args) = parser.parse_args()
    run = int(options.run)
    dir = str(options.dir)
    use_align = bool(options.al)

    tg = TransparentGrid(dir=dir, run=run)
    tg.CheckFoldersAndFiles()
    tg.OpenFileAndGetTree()
    tg.FindDiamondChannelLimits()
    tg.SetLines(try_align=use_align)
