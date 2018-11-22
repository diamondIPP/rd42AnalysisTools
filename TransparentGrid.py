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
        self.col_pitch = 50
        self.pkl = None
        self.loaded_pickle = False
        self.row_info_telescope = {'0': float(61.84669791829917), 'm': float(0.02551248435536136), 'num': 27, 'pitch': 1.008}
        self.row_info_predicted = {'0': float(61.84669791829917), 'm': float(0.02551248435536136), 'num': 27, 'pitch': 1.008}
        self.row_info_diamond = {'num': 27, 'pitch': 50}
        self.vertical_lines_telescope = []
        self.vertical_lines_telescope_tline = []
        self.horizontal_lines_telescope = []
        self.horizontal_lines_telescope_tline = []
        self.vertical_lines_diamond = []
        self.vertical_lines_diamond_tline = []
        self.horizontal_lines_diamond = []
        self.horizontal_lines_diamond_tline = []
        self.bins_per_ch_x = 3
        self.bins_per_ch_y = 3
        self.canvas = {}
        self.profile = {}
        self.histo = {}
        self.names = []
        self.tcutgs_telescope = {}
        self.tcutgs_diamond = {}
        self.goodAreas_telescope = []
        self.goodAreas_diamond = []
        self.badAreas_telescope = []
        self.badAreas_diamond = []
        self.goodAreasCutNames_telescope = ''
        self.badAreasCutNames_telescope = ''
        self.goodAreasCutNames_diamond = ''
        self.badAreasCutNames_diamond = ''

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
        self.ch_ini = int(temph.GetBinCenter(temph.FindFirstBinAbove(1, 1)))
        self.ch_end = int(temph.GetBinCenter(temph.FindLastBinAbove(1, 1)))
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
            self.row_info_telescope = self.pkl['row_info_telescope']
            self.row_info_predicted = self.pkl['row_info_telescope']
            self.vertical_lines_telescope = self.pkl['vertical_lines_telescope']
            self.horizontal_lines_telescope = self.pkl['horizontal_lines_telescope']
            self.align_info = self.pkl['align_info']
            self.row_info_diamond = self.pkl['row_info_diamond']
        elif try_align:
            self.FindHorizontalParametersThroughAlignment()
            # self.CreateLines()
        else:
            self.AskUserLowerYLines()
            self.CreateLines()

    def AskUserLowerYLines(self):
        do_telescope = raw_input('Enter 1 if you want to enter the lower y line parameters for telescope plots: ')
        if bool(int(do_telescope)):
            self.AskUserLowerYLineTelescope()
        do_diamond = raw_input('Enter 1 if you want to enter the lower y line parameters of the plots in diamond space')
        if bool(int(do_diamond)):
            self.AskUserDiamondLineParameters()

    def AskUserLowerYLineTelescope(self):
        self.row_info_telescope['0'] = self.GetFromUser('Enter the y axis intercept in silicon space for the lower detector limit (scalar between 0 and 255): ', typ='float', limits=[0, 255])
        self.row_info_telescope['m'] = self.GetFromUser('Enter the slope for the lower detector line (scalar between -1 and 1): ', typ='float', limits=[-1, 1])
        self.row_info_telescope['num'] = self.GetFromUser('Enter the number of rows: ', typ='int', limits=[1, 1000])
        self.row_info_telescope['pitch'] = self.GetFromUser('Enter the effective pitch in sil space in Y axis: ', typ='float', limits=[0, 255])

    def AskUserDiamondLineParameters(self):
        self.row_info_diamond['0'] = self.GetFromUser('Enter the y axis intercept in silicon space for the lower detector limit (scalar between 0 and 12800 in um): ', typ='float', limits=[0, 12800])
        self.row_info_diamond['pitch'] = self.GetFromUser('Enter the effective vertical pitch in um: ', typ='float', limits=[0, 20000])
        self.row_info_diamond['x_off'] = self.GetFromUser('Enter the offset for dia X ch for overlay plots (scalar between -1 and 1): ', typ='float', limits=[-1, 1])
        self.row_info_diamond['y_off'] = self.GetFromUser('Enter the offset for dia Y in um for overlay plots (scalar between -{p} and {p}): '.format(p=self.row_info_diamond['pitch']), typ='float', limits=[-self.row_info_diamond['pitch'], self.row_info_diamond['pitch']])
        self.row_info_diamond['num'] = self.GetFromUser('Enter the number of rows: ', typ='int', limits=[1, 1000])

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
        linev = self.GetVerticalLineTelescope(x=self.ch_ini - 0.5)
        lineh = self.GetHorizontalLineTelescope(y=self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini - 0.5))
        self.vertical_lines_telescope.append(linev)
        self.horizontal_lines_telescope.append(lineh)
        for col in xrange(self.num_cols):
            linev = self.GetVerticalLineTelescope(self.ch_ini + col + 0.5)
            self.vertical_lines_telescope.append(linev)
        for row in xrange(self.row_info_telescope['num']):
            lineh = self.GetHorizontalLineTelescope(y=self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini - 0.5) + (row + 1) * self.row_info_telescope['pitch'])
            self.horizontal_lines_telescope.append(lineh)

        linev = self.GetVerticalLineDiamond(x=self.ch_ini - 0.5)
        lineh = self.GetHorizontalLineDiamond(y=self.row_info_diamond['0'])
        self.vertical_lines_diamond.append(linev)
        self.horizontal_lines_diamond.append(lineh)
        for col in xrange(self.num_cols):
            linev = self.GetVerticalLineDiamond(x=self.ch_ini + col + 0.5)
            self.vertical_lines_diamond.append(linev)
        for row in xrange(self.row_info_diamond['num']):
            lineh = self.GetHorizontalLineDiamond(y=self.row_info_diamond['0'] + (row + 1) * self.row_info_diamond['pitch'])
            self.horizontal_lines_diamond.append(lineh)


    def GetVerticalLineTelescope(self, x):
        return {0: {'x': x, 'y': self.row_info_telescope['0'] + self.row_info_telescope['m'] * x}, 1: {'x': x, 'y': self.row_info_telescope['0'] + self.row_info_telescope['m'] * x + self.row_info_telescope['num'] * self.row_info_telescope['pitch']}}

    def GetHorizontalLineTelescope(self, y):
        return {0: {'x': self.ch_ini - 0.5, 'y': y}, 1: {'x': self.ch_end + 0.5, 'y': y + self.row_info_telescope['m'] * (self.ch_end - self.ch_ini + 1)}}

    def GetVerticalLineDiamond(self, x):
        return {0: {'x': x, 'y': self.row_info_diamond['0']}, 1: {'x': x, 'y': self.row_info_telescope['0'] + self.row_info_telescope['num'] * self.row_info_telescope['pitch']}}

    def GetHorizontalLineDiamond(self, y):
        return {0: {'x': self.ch_ini - 0.5, 'y': y}, 1: {'x': self.ch_end + 0.5, 'y': y}}

    def CreateLinesTLine(self):
        for lineh in self.horizontal_lines_telescope:
            self.horizontal_lines_telescope_tline.append(ro.TLine(lineh[0]['x'], lineh[0]['y'], lineh[1]['x'], lineh[1]['y']))
            self.horizontal_lines_telescope_tline[-1].SetLineColor(ro.kRed)
        for linev in self.vertical_lines_telescope:
            self.vertical_lines_telescope_tline.append(ro.TLine(linev[0]['x'], linev[0]['y'], linev[1]['x'], linev[1]['y']))
            self.vertical_lines_telescope_tline[-1].SetLineColor(ro.kRed)
        for lineh in self.horizontal_lines_diamond:
            self.horizontal_lines_diamond_tline.append(ro.TLine(lineh[0]['x'], lineh[0]['y'], lineh[1]['x'], lineh[1]['y']))
            self.horizontal_lines_diamond_tline[-1].SetLineColor(ro.kRed)
        for linev in self.vertical_lines_diamond:
            self.vertical_lines_diamond_tline.append(ro.TLine(linev[0]['x'], linev[0]['y'], linev[1]['x'], linev[1]['y']))
            self.vertical_lines_diamond_tline[-1].SetLineColor(ro.kRed)

    def DrawProfile2DFiducial(self, name, varz='clusterChargeN', cuts=''):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(128 * self.bins_per_ch_x), -0.5, 127.5, int(256 * self.bins_per_ch_y), -0.5, 255.5)
        self.profile[name].GetXaxis().SetTitle('dia X ch')
        self.profile[name].GetYaxis().SetTitle('sil ch')
        self.profile[name].GetZaxis().SetTitle('PH[ADC]')
        self.canvas[name] = ro.TCanvas('c_'+name, 'c_'+name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}:fidY:diaChXPred>>h_{n}'.format(z=varz, n=name), 'transparentEvent&&({c})'.format(c=cuts if cuts != '' else 1), 'colz prof')

    def DrawProfile2DPredicted(self, name, varz='clusterChargeN', cuts=''):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(128 * self.bins_per_ch_x), -0.5, 127.5, int(256 * self.bins_per_ch_y), 0, 12800)
        self.profile[name].GetXaxis().SetTitle('dia X ch')
        self.profile[name].GetYaxis().SetTitle('sil pred [#mum]')
        self.profile[name].GetZaxis().SetTitle('PH[ADC]')
        self.canvas[name] = ro.TCanvas('c_'+name, 'c_'+name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}:yPredicted:diaChXPred>>h_{n}'.format(z=varz, n=name), 'transparentEvent&&({c})'.format(c=cuts if cuts != '' else 1), 'colz prof')

    def DrawProfile2DDiamond(self, name, varz='clusterChargeN', cuts=''):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(128 * self.bins_per_ch_x), -0.5, 127.5, int(256 * self.bins_per_ch_y), 0, 12800)
        self.profile[name].GetXaxis().SetTitle('dia X ch')
        self.profile[name].GetYaxis().SetTitle('dia [#mum]')
        self.profile[name].GetZaxis().SetTitle('PH[ADC]')
        self.canvas[name] = ro.TCanvas('c_'+name, 'c_'+name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}:diaChYPred:diaChXPred>>h_{n}'.format(z=varz, n=name), 'transparentEvent&&({c})'.format(c=cuts if cuts != '' else 1), 'colz prof')

    def DrawLinesFiducial(self, name):
        self.DrawLines(name, type='fidY')

    def DrawLinesDiamond(self, name):
        self.DrawLines(name, type='diamond')

    def DrawLines(self, name, type):
        self.canvas[name].cd()
        if type == "fidY":
            for lineh in self.horizontal_lines_telescope_tline:
                lineh.Draw('same')
            for linev in self.vertical_lines_telescope_tline:
                linev.Draw('same')
        elif type == 'diamond':
            for lineh in self.horizontal_lines_diamond_tline:
                lineh.Draw('same')
            for linev in self.vertical_lines_diamond_tline:
                linev.Draw('same')

    def ResetLines(self):
        self.horizontal_lines_telescope = []
        self.horizontal_lines_telescope_tline = []
        self.vertical_lines_telescope = []
        self.vertical_lines_telescope_tline = []
        self.horizontal_lines_diamond = []
        self.horizontal_lines_diamond_tline = []
        self.vertical_lines_diamond = []
        self.vertical_lines_diamond_tline = []

    def CreateTCutGsTelescope(self):
        def GetNumpyArraysX(coli):
            x0 = self.ch_ini - 0.5 + coli
            x1 = self.ch_ini + 0.5 + coli
            return np.array((x0, x0, x1, x1, x0), 'f8')
        def GetNumpyArraysY(coli, rowi):
            y0 = self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini - 0.5 + coli) + rowi * self.row_info_telescope['pitch']
            y1 = self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini - 0.5 + coli) + (rowi + 1) * self.row_info_telescope['pitch']
            y2 = self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini + 0.5 + coli) + (rowi + 1) * self.row_info_telescope['pitch']
            y3 = self.row_info_telescope['0'] + self.row_info_telescope['m'] * (self.ch_ini + 0.5 + coli) + rowi * self.row_info_telescope['pitch']
            return np.array((y0, y1, y2, y3, y0), 'f8')
        for col in xrange(self.num_cols):
            self.tcutgs_telescope[col] = {}
            for row in xrange(self.row_info_telescope['num']):
                tempx = GetNumpyArraysX(col)
                tempy = GetNumpyArraysY(col, row)
                self.tcutgs_telescope[col][row] = ro.TCutG('cutg_tel_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
                self.tcutgs_telescope[col][row].SetNameTitle('cutg_tel_{c}_{r}'.format(c=col, r=row), 'cutg_tel_{c}_{r}'.format(c=col, r=row))
                self.tcutgs_telescope[col][row].SetVarX('diaChXPred')
                self.tcutgs_telescope[col][row].SetVarY('fidY')
                self.tcutgs_telescope[col][row].SetLineColor(ro.kRed)

    def CreateTCutGsDiamond(self):
        def GetNumpyArraysX(coli):
            x0 = self.ch_ini - 0.5 + coli
            x1 = self.ch_ini + 0.5 + coli
            return np.array((x0, x0, x1, x1, x0), 'f8')
        def GetNumpyArraysY(rowi):
            y0 = self.row_info_diamond['0'] + rowi * self.row_info_diamond['pitch']
            y1 = self.row_info_diamond['0'] + (rowi + 1) * self.row_info_diamond['pitch']
            return np.array((y0, y1, y1, y0, y0), 'f8')
        for col in xrange(self.num_cols):
            self.tcutgs_diamond[col] = {}
            for row in xrange(self.row_info_diamond['num']):
                tempx = GetNumpyArraysX(col)
                tempy = GetNumpyArraysY(row)
                self.tcutgs_diamond[col][row] = ro.TCutG('cutg_dia_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
                self.tcutgs_diamond[col][row].SetNameTitle('cutg_dia_{c}_{r}'.format(c=col, r=row), 'cutg_dia_{c}_{r}'.format(c=col, r=row))
                self.tcutgs_diamond[col][row].SetVarX('diaChXPred')
                self.tcutgs_diamond[col][row].SetVarY('diaChYPred')
                self.tcutgs_diamond[col][row].SetLineColor(ro.kRed)

    def AddGoodAreas(self, col, row):
        self.tcutgs_telescope[col][row].SetLineColor(ro.kRed)
        self.tcutgs_diamond[col][row].SetLineColor(ro.kRed)
        self.goodAreas_telescope.append(self.tcutgs_telescope[col][row])
        self.goodAreas_diamond.append(self.tcutgs_diamond[col][row])

    def AddBadAreas(self, col, row):
        self.tcutgs_telescope[col][row].SetLineColor(ro.kBlue)
        self.tcutgs_diamond[col][row].SetLineColor(ro.kBlue)
        self.badAreas_telescope.append(self.tcutgs_telescope[col][row])
        self.badAreas_diamond.append(self.tcutgs_diamond[col][row])

    def DrawGoodAreasTelescope(self, name):
        self.DrawGoodAreas(name, type='fidY')

    def DrawGoodAreasDiamond(self, name):
        self.DrawGoodAreas(name, type='diamond')

    def DrawGoodAreas(self, name, type):
        self.canvas[name].cd()
        if type == 'fidY':
            for area in self.goodAreas_telescope:
                area.Draw('same')
        elif type == 'diamond':
            for area in self.goodAreas_diamond:
                area.Draw('same')

    def DrawBadAreasTelescope(self, name):
        self.DrawGoodAreas(name, type='fidY')

    def DrawBadAreasDiamond(self, name):
        self.DrawGoodAreas(name, type='diamond')

    def DrawBadAreas(self, name, type):
        self.canvas[name].cd()
        if type == 'fidY':
            for area in self.badAreas_telescope:
                area.Draw('same')
        elif type == 'diamond':
            for area in self.badAreas_diamond:
                area.Draw('same')

    def SelectGoodAndBadByThreshold(self, val=500):
        for col in xrange(self.num_cols):
            for row in xrange(self.row_info_diamond['num']):
                temph = ro.TH1F('temphrc', 'temphrc', 200, 0, 4000)
                self.trans_tree.Draw('clusterChargeN>>temphrc', 'transparentEvent&&({n})'.format(n=self.tcutgs_diamond[col][row].GetName()), 'goff')
                if temph.GetMean() > val:
                    self.AddGoodAreas(col, row)
                else:
                    self.AddBadAreas(col, row)
                del temph

        tempgood = [cut.GetName() for cut in self.goodAreas_telescope]
        self.goodAreasCutNames_telescope = '((' + ')||('.join(tempgood) + '))'
        tempbad = [cut.GetName() for cut in self.badAreas_telescope]
        self.badAreasCutNames_telescope = '((' + ')||('.join(tempbad) + '))'
        
        tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
        self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
        tempbad = [cut.GetName() for cut in self.badAreas_diamond]
        self.badAreasCutNames_diamond = '((' + ')||('.join(tempbad) + '))'

    def ResetAreas(self):
        self.goodAreas_telescope = []
        self.badAreas_telescope = []
        self.goodAreas_diamond = []
        self.badAreas_diamond = []

    def DrawPHGoodAreas(self, name, var='clusterChargeN', type='diamond', cuts=''):
        self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, self.phbins, self.phmin, self.phmax)
        self.histo[name].GetXaxis().SetTitle('PH[ADC]')
        self.histo[name].GetYaxis().SetTitle('entries')
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}>>h_{n}'.format(z=var, n=name), 'transparentEvent&&{n}&&({c})'.format(n=self.goodAreasCutNames_diamond if type=='diamond' else self.goodAreasCutNames_telescope, c=cuts if cuts != '' else 1))

    def DrawPHBadAreas(self, name, var='clusterChargeN', type='diamond', cuts=''):
        self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, self.phbins, self.phmin, self.phmax)
        self.histo[name].GetXaxis().SetTitle('PH[ADC]')
        self.histo[name].GetYaxis().SetTitle('entries')
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{z}>>h_{n}'.format(z=var, n=name), 'transparentEvent&&{n}&&({c})'.format(n=self.badAreasCutNames_diamond if type=='diamond' else self.badAreasCutNames_telescope, c=cuts if cuts != '' else 1))

    def Draw2DProfileDiamondChannelOverlay(self, name, var='clusterChargeN', cells='all', cut=''):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(self.col_pitch + 10), -5, self.col_pitch + 5, int(256 * self.bins_per_ch_y), 0, 12800)
        self.profile[name].GetXaxis().SetTitle('dia X [#mum]')
        self.profile[name].GetYaxis().SetTitle('dia Y [#mum]')
        self.profile[name].GetZaxis().SetTitle('PH[ADC]')
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        temp_cuts = '{n}'.format(n=self.goodAreasCutNames_diamond) if cells == 'good' else '{n}'.format(n=self.badAreasCutNames_diamond) if cells == 'bad' else '(1)'
        temp_cuts = temp_cuts if cut == '' else temp_cuts + '&&({c})'.format(c=cut)
        self.trans_tree.Draw('{z}:diaChYPred:((diaChXPred-{o})*{p})%{p}>>h_{n}'.format(z=var, o=self.row_info_diamond['x_off'], n=name, p=self.col_pitch), 'transparentEvent&&{c}'.format(c=temp_cuts), 'colz prof')

    def Draw2DProfileDiamondCellOverlay(self, name, var='clusterChargeN', cells='all', cut=''):
        self.profile[name] = ro.TProfile2D('h_'+name, 'h_'+name, int(self.col_pitch + 10), -5, self.col_pitch + 5, int(self.row_info_diamond['pitch']) + 10, -5, int(self.row_info_diamond['pitch']) + 5)
        self.profile[name].GetXaxis().SetTitle('dia X [#mum]')
        self.profile[name].GetYaxis().SetTitle('dia Y [#mum]')
        self.profile[name].GetZaxis().SetTitle('PH[ADC]')
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        temp_cuts = '({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['0'] + self.row_info_diamond['pitch'] * self.row_info_diamond['num'])
        temp_cuts += '{n}'.format(n=self.goodAreasCutNames_diamond) if cells == 'good' else '{n}'.format(n=self.badAreasCutNames_diamond) if cells == 'bad' else '(1)'
        temp_cuts = temp_cuts if cut == '' else temp_cuts + '&&({c})'.format(c=cut)
        self.trans_tree.Draw('{z}:(((diaChYPred-{oy})*10000)%{srp})/10000:((diaChXPred-{ox})*{p})%{p}>>h_{n}'.format(z=var, ox=self.row_info_diamond['x_off'], n=name, oy=self.row_info_diamond['y_off'], srp=int(10000*self.row_info_diamond['pitch']), p=self.col_pitch), 'transparentEvent&&{c}'.format(c=temp_cuts), 'colz prof')


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
