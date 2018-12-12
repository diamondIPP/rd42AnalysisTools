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
from Langaus import LanGaus

class TransparentGrid:
    def __init__(self, dir='', run=25209, cellsize=50.0):
        ro.gStyle.SetPalette(55)
        ro.gStyle.SetNumberContours(99)
        ro.gStyle.SetOptStat(11)
        ro.gStyle.SetOptFit(111)
        ro.gStyle.SetPaintTextFormat(".0f")
        # ro.TFormula.SetMaxima(100000)
        self.run = run
        self.dir = os.path.abspath(os.path.expanduser(os.path.expandvars(dir)))
        self.trans_file = None
        self.trans_tree = None
        self.align_file = None
        self.align_obj = None
        self.align_info = {'xoff': float(0), 'phi': float(0)}
        self.num_cols = 19 if cellsize == 50 else 13
        self.ch_ini = 0
        self.ch_end = 84
        self.phbins = 200
        self.phmin = 0
        self.phmax = 4000
        self.neg_cut = 3
        self.col_pitch = cellsize
        self.cell_resolution = 50.0 / 25 if self.col_pitch == 50 else 100.0 / 51
        self.saturated_ADC = 4095
        self.pkl = None
        self.loaded_pickle = False
        self.row_info_telescope = {'0': float(61.95469791829917), 'm': float(0.02551248435536136), 'num': 27, 'pitch': 1.0} if self.col_pitch == 50 else {'0': float(61.95469791829917), 'm': float(0.05102496), 'num': 24, 'pitch': 1.0}
        self.row_info_predicted = {'0': float(61.95469791829917), 'm': float(0.02551248435536136), 'num': 27, 'pitch': 1.0} if self.col_pitch == 50 else {'0': float(61.95469791829917), 'm': float(0.05102496), 'num': 24, 'pitch': 1.0}
        if self.run in [25207, 25208, 25209, 25210]:
            if self.run in [25207, 25208, 25210]: print 'Using settings for run 25209. Results might not be correctly aligned'
            self.row_info_diamond = {'num': 27, 'pitch': 50.0, 'x_off': 0.5065, 'y_off': 18.95, '0': 3117.70005, 'up': 4467.70005} if self.col_pitch == 50 else {'num': 24, 'pitch': 100.0, 'x_off': 0.4976, 'y_off': 52.6, '0': 3049.6, 'up': 5449.6}
        elif self.run in [25202, 25203, 25204, 25205, 25206]:
            if self.run in [25202, 25203, 25206]: print 'Using settings for run 25205. Results might not be correctly aligned'
            self.row_info_diamond = {'num': 27, 'pitch': 50.0, 'x_off': 0.4922, 'y_off': 47.978033, '0': 3148.785, 'up': 4498.785} if self.col_pitch == 50 else {'num': 24, 'pitch': 100.0, 'x_off': 0.4976, 'y_off': 79.803, '0': 3076.115, 'up': 5476.115}
        self.vertical_lines_telescope = []
        self.vertical_lines_telescope_tline = []
        self.horizontal_lines_telescope = []
        self.horizontal_lines_telescope_tline = []
        self.vertical_lines_diamond = []
        self.vertical_lines_diamond_tline = []
        self.horizontal_lines_diamond = []
        self.horizontal_lines_diamond_tline = []
        self.bins_per_ch_x = 3 if self.col_pitch == 50 else 5
        self.bins_per_ch_y = 3 if self.col_pitch == 50 else 5
        self.canvas = {}
        self.line = {}
        self.profile = {}
        self.histo = {}
        self.names = []
        self.tcutgs_telescope = {}
        self.tcutgs_diamond = {}
        self.tcutg_diamond_center = None
        self.tcutgs_diamond_center = {}
        self.length_central_region = 30 if self.col_pitch == 50 else 40
        self.goodAreas_telescope = []
        self.goodAreas_diamond = []
        self.goodAreas_diamond_centers = []
        self.goodAreas_index = []
        self.badAreas_telescope = []
        self.badAreas_diamond = []
        self.badAreas_diamond_centers = []
        self.badAreas_index = []
        self.goodAreasCutNames_telescope = ''
        self.badAreasCutNames_telescope = ''
        self.goodAreasCutNames_diamond = ''
        self.goodAreasCutNames_diamond_centers = ''
        self.badAreasCutNames_diamond = ''
        self.badAreasCutNames_diamond_centers = ''
        self.temph = None
        self.langaus = {}
        self.gridTextDiamond = None
        self.gridTextTelescope = None
        self.doubleLangaus = {}
        self.lg1, self.lg2 = ro.TF1(), ro.TF1()

        self.conv_steps = 1000
        self.sigma_conv = 5
        self.mpshift = -0.22278298


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
        self.trans_tree.Draw('diaChXPred>>temph', 'transparentEvent', 'goff')
        self.ch_ini = int(temph.GetBinCenter(temph.FindFirstBinAbove(1, 1)))
        self.ch_end = int(temph.GetBinCenter(temph.FindLastBinAbove(1, 1)))
        self.num_cols = self.ch_end - self.ch_ini + 1

    def TryLoadPickle(self):
        picklepath = '{d}/{r}/transp_grid.{r}.pkl'.format(d=self.dir, r=self.run)
        if os.path.isfile(picklepath):
            with open(picklepath, 'rb') as pkl:
                self.pkl = pickle.load(pkl)
                self.loaded_pickle = True

    def SetLines(self, try_align=False):
        self.TryLoadPickle()
        if self.loaded_pickle:
            self.row_info_telescope = self.pkl['row_info_telescope']
            self.row_info_predicted = self.pkl['row_info_telescope']
            self.row_info_diamond = self.pkl['row_info_diamond']
            self.vertical_lines_telescope = self.pkl['vertical_lines_telescope']
            self.vertical_lines_diamond = self.pkl['vertical_lines_diamond']
            self.horizontal_lines_diamond = self.pkl['horizontal_lines_diamond']
            self.num_cols = self.pkl['num_cols']
            self.col_pitch = self.pkl['col_pitch']
            self.horizontal_lines_telescope = self.pkl['horizontal_lines_telescope']
            self.align_info = self.pkl['align_info']
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
        return {0: {'x': x, 'y': self.row_info_diamond['0']}, 1: {'x': x, 'y': self.row_info_diamond['0'] + self.row_info_diamond['num'] * self.row_info_diamond['pitch']}}

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

    def DrawProfile2D(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, varz='clusterChargeN', zname='PH[ADC]', cuts=''):
        # ro.gStyle.SetOptStat('en')
        ro.TFormula.SetMaxima(100000)
        self.profile[name] = ro.TProfile2D('h_' + name, 'h_' + name, int(np.floor((xmax - xmin)/deltax + 0.5) + 2), xmin - deltax, xmax + deltax, int(np.floor((ymax - ymin)/deltay + 0.5) + 2), ymin - deltay, ymax + deltay)
        self.profile[name].GetXaxis().SetTitle(xname)
        self.profile[name].GetYaxis().SetTitle(yname)
        self.profile[name].GetZaxis().SetTitle(zname)
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        list_cuts = ['transparentEvent']
        if cuts != '':
            list_cuts.append(cuts)
        temp_cut = '&&'.join(list_cuts)
        self.trans_tree.Draw('{z}:{y}:{x}>>h_{n}'.format(z=varz, y=vary, x=varx, n=name), temp_cut, 'colz prof')
        ro.gPad.Update()
        if self.profile[name].FindObject('stats'):
            self.profile[name].FindObject('stats').SetOptStat(11)
            self.profile[name].FindObject('stats').SetX1NDC(0.7)
            self.profile[name].FindObject('stats').SetX2NDC(0.9)
            self.profile[name].FindObject('stats').SetY1NDC(0.9)
            self.profile[name].FindObject('stats').SetY2NDC(0.975)
        if self.profile[name].FindObject('palette'):
            self.profile[name].FindObject('palette').SetX1NDC(0.87)
            self.profile[name].FindObject('palette').SetX2NDC(0.92)
        ro.gPad.Update()
        ro.TFormula.SetMaxima(1000)

    def Draw2DHisto(self, name, xmin, xmax, deltax, xname, ymin, ymax, deltay, yname, varx, vary, cuts=''):
        # ro.gStyle.SetOptStat('en')
        ro.TFormula.SetMaxima(100000)
        self.histo[name] = ro.TH2F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5) + 2), xmin - deltax, xmax + deltax, int(np.floor((ymax - ymin) / deltay + 0.5) + 2), ymin - deltay, ymax + deltay)
        self.histo[name].GetXaxis().SetTitle(xname)
        self.histo[name].GetYaxis().SetTitle(yname)
        self.histo[name].GetZaxis().SetTitle('entries')
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        list_cuts = ['transparentEvent']
        if cuts != '':
            list_cuts.append(cuts)
        temp_cut = '&&'.join(list_cuts)
        self.trans_tree.Draw('{y}:{x}>>h_{n}'.format(y=vary, x=varx, n=name), temp_cut, 'colz')
        ro.gPad.Update()
        if self.histo[name].FindObject('stats'):
            self.histo[name].FindObject('stats').SetOptStat(11)
            self.histo[name].FindObject('stats').SetX1NDC(0.7)
            self.histo[name].FindObject('stats').SetX2NDC(0.9)
            self.histo[name].FindObject('stats').SetY1NDC(0.9)
            self.histo[name].FindObject('stats').SetY2NDC(0.975)
        if self.histo[name].FindObject('palette'):
            self.histo[name].FindObject('palette').SetX1NDC(0.87)
            self.histo[name].FindObject('palette').SetX2NDC(0.92)
        ro.gPad.Update()
        ro.TFormula.SetMaxima(1000)

    def DrawPH(self, name, xmin, xmax, deltax, var='clusterChargeN', varname='PH[ADC]', cuts=''):
        ro.TFormula.SetMaxima(100000)
        # ro.gStyle.SetOptStat('neMmRruo')
        self.histo[name] = ro.TH1F('h_' + name, 'h_' + name, int(np.floor((xmax - xmin) / deltax + 0.5)), xmin, xmax)
        self.histo[name].GetXaxis().SetTitle(varname)
        self.histo[name].GetYaxis().SetTitle('entries')
        list_cuts = ['transparentEvent']
        if cuts != '':
            list_cuts.append(cuts)
        temp_cuts = '&&'.join(list_cuts)
        self.canvas[name] = ro.TCanvas('c_' + name, 'c_' + name, 1)
        self.canvas[name].cd()
        self.trans_tree.Draw('{v}>>h_{n}'.format(v=var, n=name), temp_cuts, 'e')
        self.canvas[name].SetGridx()
        self.canvas[name].SetGridy()
        self.canvas[name].SetTicky()
        ro.gPad.Update()
        self.histo[name].FindObject('stats').SetOptStat(112211)
        self.histo[name].FindObject('stats').SetX1NDC(0.6)
        self.histo[name].FindObject('stats').SetX2NDC(0.9)
        self.histo[name].FindObject('stats').SetY1NDC(0.6)
        self.histo[name].FindObject('stats').SetY2NDC(0.9)
        ro.gPad.Update()
        ro.TFormula.SetMaxima(1000)

    def DrawProfile2DFiducial(self, name, varz='clusterChargeN', cuts=''):
        self.DrawProfile2D(name, -0.5, 127.5, 1.0/self.bins_per_ch_x, 'dia X ch', -0.5, 255.5, float(self.row_info_telescope['pitch'])/self.bins_per_ch_y, 'sil Y ch', 'diaChXPred', 'fidY', varz, 'PH[ADC]', cuts)

    def DrawProfile2DPredicted(self, name, varz='clusterChargeN', cuts=''):
        self.DrawProfile2D(name, -0.5, 127.5, 1.0/self.bins_per_ch_x, 'dia X ch', 0, 12800, 50.0 * self.row_info_predicted['pitch']/self.bins_per_ch_y, 'sil pred Y [#mum]', 'diaChXPred', 'yPredicted', varz, 'PH[ADC]', cuts)

    def DrawProfile2DDiamond(self, name, varz='clusterChargeN', cuts='', draw_top_borders=False):
        list_cuts = []
        if not draw_top_borders:
            list_cuts.append('({l}<=diaChYPred)&&(diaChYPred<={h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['up']))
        if cuts != '':
            list_cuts.append(cuts)
        temp_cuts = '&&'.join(list_cuts)
        self.DrawProfile2D(name, -0.5, 127.5, 1.0/self.bins_per_ch_x, 'dia X ch', self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'], float(self.row_info_diamond['pitch'])/self.bins_per_ch_y, 'sil pred Y [#mum]', 'diaChXPred', 'diaChYPred', varz, 'PH[ADC]', temp_cuts)

    def GetOccupancyFromProfile(self, name):
        # ro.gStyle.SetOptStat('ne')
        name_occupancy = 'hit_map_' + name
        self.histo[name_occupancy] = self.profile[name].ProjectionXY('h_' + name_occupancy, 'B')
        self.histo[name_occupancy].SetTitle('h_' + name_occupancy)
        self.histo[name_occupancy].GetXaxis().SetTitle(self.profile[name].GetXaxis().GetTitle())
        self.histo[name_occupancy].GetYaxis().SetTitle(self.profile[name].GetYaxis().GetTitle())
        self.histo[name_occupancy].GetZaxis().SetTitle('entries')
        self.canvas[name_occupancy] = ro.TCanvas('c_' + name_occupancy, 'c_' + name_occupancy, 1)
        self.canvas[name_occupancy].cd()
        self.histo[name_occupancy].Draw('colz')
        ro.gPad.Update()
        self.histo[name_occupancy].FindObject('stats').SetOptStat(11)
        self.histo[name_occupancy].FindObject('stats').SetX1NDC(0.7)
        self.histo[name_occupancy].FindObject('stats').SetX2NDC(0.9)
        self.histo[name_occupancy].FindObject('stats').SetY1NDC(0.9)
        self.histo[name_occupancy].FindObject('stats').SetY2NDC(0.975)
        ro.gPad.Update()

    def DrawLinesFiducial(self, name):
        self.DrawLines(name, type='fidY')

    def DrawLinesDiamond(self, name):
        self.DrawLines(name, type='diamond')

    def DrawLines(self, name, type):
        # ro.gStyle.SetOptStat('en')
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

    def CreateTCutGs(self):
        self.CreateTCutGsTelescope()
        self.CreateTCutGsDiamond()
        self.CreateTCutGsDiamondCenter()
        self.CreateGridText()

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
                self.tcutgs_telescope[col][row].SetLineColor(ro.kBlack)

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
                self.tcutgs_diamond[col][row].SetLineColor(ro.kBlack)
                if col == self.num_cols -1 and row == self.row_info_diamond['num'] - 1:
                    self.row_info_diamond['up'] = tempy[2]

    def CreateTCutGsDiamondCenter(self):
        def GetNumpyArraysX(coli):
            x0 = self.ch_ini - self.length_central_region/(2.0*self.col_pitch) + coli
            x1 = self.ch_ini + self.length_central_region/(2.0*self.col_pitch) + coli
            return np.array((x0, x0, x1, x1, x0), 'f8')

        def GetNumpyArraysY(rowi):
            y0 = self.row_info_diamond['0'] + rowi * self.row_info_diamond['pitch'] + self.row_info_diamond['pitch']/2.0 - self.length_central_region/2.0
            y1 = self.row_info_diamond['0'] + (rowi + 1) * self.row_info_diamond['pitch'] - self.row_info_diamond['pitch']/2.0 + self.length_central_region/2.0
            return np.array((y0, y1, y1, y0, y0), 'f8')

        x0i = self.col_pitch / 2.0 - self.length_central_region / 2.0
        x1i = self.col_pitch / 2.0 + self.length_central_region / 2.0
        y0i = self.row_info_diamond['pitch'] / 2.0 - self.length_central_region / 2.0
        y1i = self.row_info_diamond['pitch'] / 2.0 + self.length_central_region / 2.0
        tempi = np.array((x0i, x0i, x1i, x1i, x0i), 'f8')
        tempj = np.array((y0i, y1i, y1i, y0i, y0i), 'f8')
        self.tcutg_diamond_center = ro.TCutG('cutg_dia_center', 5, tempi, tempj)
        self.tcutg_diamond_center.SetNameTitle('cutg_dia_center', 'cutg_dia_center')
        self.tcutg_diamond_center.SetVarX('diaChXPred')
        self.tcutg_diamond_center.SetVarY('diaChYPred')
        self.tcutg_diamond_center.SetLineColor(ro.kViolet)

        for col in xrange(self.num_cols):
            self.tcutgs_diamond_center[col] = {}
            for row in xrange(self.row_info_diamond['num']):
                tempx = GetNumpyArraysX(col)
                tempy = GetNumpyArraysY(row)
                self.tcutgs_diamond_center[col][row] = ro.TCutG('cutg_dia_center_{c}_{r}'.format(c=col, r=row), 5, tempx, tempy)
                self.tcutgs_diamond_center[col][row].SetNameTitle('cutg_dia_center_{c}_{r}'.format(c=col, r=row), 'cutg_center_dia_{c}_{r}'.format(c=col, r=row))
                self.tcutgs_diamond_center[col][row].SetVarX('diaChXPred')
                self.tcutgs_diamond_center[col][row].SetVarY('diaChYPred')
                self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)

    def CreateGridText(self):
        self.gridTextDiamond = ro.TH2F('gridText_diamond', 'gridText_diamond', int(np.floor(128.0 * self.bins_per_ch_x + 0.5) + 2), -0.5 - 1.0 / self.bins_per_ch_x, 127.5 + 1.0 / self.bins_per_ch_x, int(np.floor(256 * self.bins_per_ch_y + 0.5) + 2), self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'] - (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y), self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'] + (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y))
        x0, x1, y0, y1 = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
        for col in xrange(0, self.num_cols):
            self.tcutgs_diamond[col][0].GetPoint(0, x0, y0)
            self.tcutgs_diamond[col][0].GetPoint(3, x1, y0)
            self.gridTextDiamond.Fill(np.mean((x0, x1)), y0[0]-0.1, (col + 0.01))
        for row in xrange(0, self.row_info_diamond['num']):
            self.tcutgs_diamond[0][row].GetPoint(0, x0, y0)
            self.tcutgs_diamond[0][row].GetPoint(1, x0, y1)
            self.gridTextDiamond.Fill(x0[0]-0.1, np.mean((y0, y1)), (row + 0.01))
        self.gridTextDiamond.SetMarkerSize(0.8)
        self.gridTextTelescope = ro.TH2F('gridText_telescope', 'gridText_telescope', int(np.floor(128.0 * self.bins_per_ch_x + 0.5) + 2), -0.5 - 1.0 / self.bins_per_ch_x, 127.5 + 1.0 / self.bins_per_ch_x, int(np.floor(256 * self.bins_per_ch_y + 0.5) + 2), self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'] - (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y), self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'] + (float(self.row_info_diamond['pitch']) / self.bins_per_ch_y))
        x0, x1, y0, y1 = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
        for col in xrange(0, self.num_cols):
            self.tcutgs_telescope[col][0].GetPoint(0, x0, y0)
            self.tcutgs_telescope[col][0].GetPoint(3, x1, y0)
            self.gridTextTelescope.Fill(np.mean((x0, x1)), y0[0]-0.1, (col + 0.01))
        for row in xrange(0, self.row_info_telescope['num']):
            self.tcutgs_telescope[0][row].GetPoint(0, x0, y0)
            self.tcutgs_telescope[0][row].GetPoint(1, x0, y1)
            self.gridTextTelescope.Fill(x0[0]-0.1, np.mean((y0, y1)), (row + 0.01))
        self.gridTextTelescope.SetMarkerSize(0.8)

    def AddGoodAreas(self, col, row):
        if 0 <= col < self.num_cols and 0 <= row < self.row_info_diamond['num']:
            self.goodAreas_index.append((col, row))
            self.tcutgs_telescope[col][row].SetLineColor(ro.kRed)
            self.tcutgs_diamond[col][row].SetLineColor(ro.kRed)
            self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
            self.goodAreas_telescope.append(self.tcutgs_telescope[col][row])
            self.goodAreas_diamond.append(self.tcutgs_diamond[col][row])
            self.goodAreas_diamond_centers.append(self.tcutgs_diamond_center[col][row])

            tempgood = [cut.GetName() for cut in self.goodAreas_telescope]
            self.goodAreasCutNames_telescope = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
            self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
            self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

    def AddGoodAreasRow(self, row, coli=0, colf=0):
        (colii, colff) = (0, self.num_cols) if coli == 0 and colf == 0 else (coli, colf)
        if 0 <= colii <= colff < self.num_cols and 0 <= row < self.row_info_diamond['num']:
            for col in xrange(colii, colff + 1):
                self.goodAreas_index.append((col, row))
                self.tcutgs_telescope[col][row].SetLineColor(ro.kRed)
                self.tcutgs_diamond[col][row].SetLineColor(ro.kRed)
                self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
                self.goodAreas_telescope.append(self.tcutgs_telescope[col][row])
                self.goodAreas_diamond.append(self.tcutgs_diamond[col][row])
                self.goodAreas_diamond_centers.append(self.tcutgs_diamond_center[col][row])

            tempgood = [cut.GetName() for cut in self.goodAreas_telescope]
            self.goodAreasCutNames_telescope = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
            self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
            self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

    def AddGoodAreasCol(self, col, rowi=0, rowf=0):
        (rowii, rowff) = (0, self.row_info_diamond['num']) if rowi == 0 and rowf == 0 else (rowi, rowf)
        if 0 <= col < self.num_cols and 0 <= rowii <= rowff < self.row_info_diamond['num']:
            for row in xrange(rowii, rowff + 1):
                self.goodAreas_index.append((col, row))
                self.tcutgs_telescope[col][row].SetLineColor(ro.kRed)
                self.tcutgs_diamond[col][row].SetLineColor(ro.kRed)
                self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
                self.goodAreas_telescope.append(self.tcutgs_telescope[col][row])
                self.goodAreas_diamond.append(self.tcutgs_diamond[col][row])
                self.goodAreas_diamond_centers.append(self.tcutgs_diamond_center[col][row])

            tempgood = [cut.GetName() for cut in self.goodAreas_telescope]
            self.goodAreasCutNames_telescope = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
            self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
            self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

    def AddBadAreas(self, col, row):
        if 0 <= col < self.num_cols and 0 <= row < self.row_info_diamond['num']:
            self.badAreas_index.append((col, row))
            self.tcutgs_telescope[col][row].SetLineColor(ro.kBlue)
            self.tcutgs_diamond[col][row].SetLineColor(ro.kBlue)
            self.tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
            self.badAreas_telescope.append(self.tcutgs_telescope[col][row])
            self.badAreas_diamond.append(self.tcutgs_diamond[col][row])
            self.badAreas_diamond_centers.append(self.tcutgs_diamond_center[col][row])

            tempbad = [cut.GetName() for cut in self.badAreas_telescope]
            self.badAreasCutNames_telescope = '((' + ')||('.join(tempbad) + '))'
            tempbad = [cut.GetName() for cut in self.badAreas_diamond]
            self.badAreasCutNames_diamond = '((' + ')||('.join(tempbad) + '))'
            tempbad = [cut.GetName() for cut in self.badAreas_diamond_centers]
            self.badAreasCutNames_diamond_centers = '((' + ')||('.join(tempbad) + '))'

    def AddRemainingToBadAreas(self):
        for col in xrange(0, self.num_cols):
            for row in xrange(0, self.row_info_diamond['num']):
                if (col, row) not in self.goodAreas_index and (col, row) not in self.badAreas_index:
                    self.AddBadAreas(col, row)

    def RemoveFromGoodArea(self, col, row):
        if (col, row) in self.goodAreas_index:
            index_g = self.goodAreas_index.index((col, row))
            self.goodAreas_telescope.pop(index_g)
            self.goodAreas_diamond.pop(index_g)
            self.goodAreas_diamond_centers.pop(index_g)
            self.goodAreas_index.pop(index_g)
            self.AddBadAreas(col, row)
            tempgood = [cut.GetName() for cut in self.goodAreas_telescope]
            self.goodAreasCutNames_telescope = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
            self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
            tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
            self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

    def DrawGoodAreasTelescope(self, name):
        self.DrawGoodAreas(name, type='fidY')

    def DrawGoodAreasDiamond(self, name):
        self.DrawGoodAreas(name, type='diamond')

    def DrawGoodAreasDiamondCenters(self, name):
        self.DrawGoodAreas(name, type='centers')

    def DrawTCutGs(self, name, type):
        self.canvas[name].cd()
        # ro.gStyle.SetOptStat('en')
        ro.gStyle.SetPaintTextFormat(".0f")
        if type == 'fidY':
            self.gridTextTelescope.Draw('same TEXT0')
        elif type == 'diamond':
            self.gridTextDiamond.Draw('same TEXT0')
        if name in self.profile.keys():
            self.profile[name].Draw('same colz')
        elif name in self.histo.keys():
            self.histo[name].Draw('same colz')
        for col in xrange(0, self.num_cols):
            for row in xrange(0, self.row_info_diamond['num']):
                if type == 'fidY':
                    self.tcutgs_telescope[col][row].Draw('same')
                elif type == 'diamond':
                    self.tcutgs_diamond[col][row].Draw('same')
                elif type == 'centers':
                    self.tcutgs_diamond_center[col][row].Draw('same')

    def DrawGoodAreas(self, name, type):
        # ro.gStyle.SetOptStat('en')
        self.canvas[name].cd()
        if type == 'fidY':
            for area in self.goodAreas_telescope:
                area.Draw('same')
        elif type == 'diamond':
            for area in self.goodAreas_diamond:
                area.Draw('same')
        elif type == 'centers':
            for area in self.goodAreas_diamond_centers:
                area.Draw('same')

    def DrawBadAreasTelescope(self, name):
        self.DrawBadAreas(name, type='fidY')

    def DrawBadAreasDiamond(self, name):
        self.DrawBadAreas(name, type='diamond')

    def DrawBadAreasDiamondCenters(self, name):
        self.DrawBadAreas(name, type='centers')

    def DrawBadAreas(self, name, type):
        # ro.gStyle.SetOptStat('en')
        self.canvas[name].cd()
        if type == 'fidY':
            for area in self.badAreas_telescope:
                area.Draw('same')
        elif type == 'diamond':
            for area in self.badAreas_diamond:
                area.Draw('same')
        elif type == 'centers':
            for area in self.badAreas_diamond_centers:
                area.Draw('same')

    def SelectGoodAndBadByThreshold(self, val=500, var='clusterChargeN'):
        for col in xrange(self.num_cols):
            for row in xrange(self.row_info_diamond['num']):
                # self.temph = ro.TH1F('temphrc', 'temphrc', 200, 0, 4000)
                self.trans_tree.Draw(var+'>>temphrc(200,0,4000)', 'transparentEvent&&({n})'.format(n=self.tcutgs_diamond[col][row].GetName()), 'goff')
                temph = ro.gDirectory.Get('temphrc')
                if temph.GetMean() > val:
                    self.AddGoodAreas(col, row)
                else:
                    self.AddBadAreas(col, row)
                temph.Reset('ICES')
                temph.Delete()
                del temph

    def ResetAreas(self):
        self.goodAreas_telescope = []
        self.badAreas_telescope = []
        self.goodAreas_diamond = []
        self.goodAreas_diamond_centers = []
        self.badAreas_diamond = []
        self.badAreas_diamond_centers = []
        self.goodAreasCutNames_diamond = ''
        self.goodAreasCutNames_diamond_centers = ''
        self.badAreasCutNames_diamond = ''
        self.badAreasCutNames_diamond_centers = ''

    def DrawPHGoodAreas(self, name, var='clusterChargeN', cuts='', type='diamond'):
        list_cuts = ['{n}'.format(n=self.goodAreasCutNames_diamond if type == 'diamond' else self.goodAreasCutNames_telescope)]
        if cuts != '':
            list_cuts.append(cuts)
        temp_cut = '&&'.join(list_cuts)
        self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cut)

    def DrawPHBadAreas(self, name, var='clusterChargeN', cuts='', type='diamond'):
        list_cuts = ['{n}'.format(n=self.badAreasCutNames_diamond if type == 'diamond' else self.badAreasCutNames_telescope)]
        if cuts != '':
            list_cuts.append(cuts)
        temp_cut = '&&'.join(list_cuts)
        self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cut)

    def DrawPHCentralRegion(self, name, var='clusterChargeN', cells='good', cuts=''):
        list_cuts = ['{n}'.format(n=self.goodAreasCutNames_diamond_centers) if cells == 'good' else '{n}'.format(n=self.badAreasCutNames_diamond_centers) if cells == 'bad' else '({n}||{m})'.format(n=self.goodAreasCutNames_diamond_centers, m=self.badAreasCutNames_diamond_centers)]
        if cuts != '':
            list_cuts.append(cuts)
        temp_cuts = '&&'.join(list_cuts)
        self.DrawPH(name, self.phmin, self.phmax, (self.phmax - self.phmin) / self.phbins, var, 'PH[ADC]', temp_cuts)

    def Draw2DProfileDiamondChannelOverlay(self, name, var='clusterChargeN', cells='all', cut=''):
        list_cuts = ['{n}'.format(n=self.goodAreasCutNames_diamond) if cells == 'good' else '{n}'.format(n=self.badAreasCutNames_diamond) if cells == 'bad' else '(1)']
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.DrawProfile2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'],
                           float(self.row_info_diamond['pitch'])/self.bins_per_ch_y, 'dia Y [#mum]', '((diaChXPred-{o})*{p})%{p}'.format(o=self.row_info_diamond['x_off'], p=self.col_pitch), 'diaChYPred', var, 'PH[ADC]', temp_cuts)

    def Draw2DProfileDiamondRowOverlay(self, name, var='clusterChargeN', cells='all', cut=''):
        list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['0'] + self.row_info_diamond['pitch'] * self.row_info_diamond['num'])]
        if cells == 'good':
            list_cuts.append(self.goodAreasCutNames_diamond)
        elif cells == 'bad':
            list_cuts.append(self.badAreasCutNames_diamond)
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.DrawProfile2D(name, -0.5, 127.5, self.cell_resolution, 'dia X ch', 0, self.row_info_diamond['pitch'], self.cell_resolution, 'dia Y [#mum]', 'diaChXPred', '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=self.row_info_diamond['y_off'], srp=int(10000*self.row_info_diamond['pitch'])), var, 'PH[ADC]', temp_cuts)

    def Draw2DProfileDiamondCellOverlay(self, name, var='clusterChargeN', cells='all', cut=''):
        list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['0'] + self.row_info_diamond['pitch'] * self.row_info_diamond['num'])]
        if cells == 'good':
            list_cuts.append(self.goodAreasCutNames_diamond)
        elif cells == 'bad':
            list_cuts.append(self.badAreasCutNames_diamond)
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.DrawProfile2D(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', 0, self.row_info_diamond['pitch'], self.cell_resolution, 'dia Y [#mum]', '((diaChXPred-{ox})*{p})%{p}'.format(ox=self.row_info_diamond['x_off'], p=self.col_pitch), '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=self.row_info_diamond['y_off'], srp=int(10000*self.row_info_diamond['pitch'])), var, 'PH[ADC]', temp_cuts)

    def Draw2DHistoDiamond(self, name, cut=''):
        self.Draw2DHisto(name, -0.5, 127.5, 1.0 / (self.bins_per_ch_x), 'dia X ch', self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'],
                         float(self.row_info_diamond['pitch']) / (self.bins_per_ch_y), 'dia Y [#mum]', 'diaChXPred', 'diaChYPred', cut)

    def Draw2DHistoDiamondChannelOverlay(self, name, cells='all', cut=''):
        list_cuts = []
        if cells == 'good':
            list_cuts.append(self.goodAreasCutNames_diamond)
        elif cells == 'bad':
            list_cuts.append(self.badAreasCutNames_diamond)
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.Draw2DHisto(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', self.row_info_diamond['0'] - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5) * self.row_info_diamond['pitch'], self.row_info_diamond['0'] + (256 - np.floor(self.row_info_diamond['0'] / self.row_info_diamond['pitch'] + 0.5)) * self.row_info_diamond['pitch'],
                         float(self.row_info_diamond['pitch']) / self.bins_per_ch_y, 'dia Y [#mum]', '((diaChXPred-{o})*{p})%{p}'.format(o=self.row_info_diamond['x_off'], p=self.col_pitch), 'diaChYPred', temp_cuts)

    def Draw2DHistoDiamondRowOverlay(self, name, cells='all', cut=''):
        list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['0'] + self.row_info_diamond['pitch'] * self.row_info_diamond['num'])]
        if cells == 'good':
            list_cuts.append(self.goodAreasCutNames_diamond)
        elif cells == 'bad':
            list_cuts.append(self.badAreasCutNames_diamond)
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.Draw2DHisto(name, -0.5, 127.5, 1.0 / self.bins_per_ch_x, 'dia X ch', 0, self.row_info_diamond['pitch'], self.cell_resolution, 'dia Y [#mum]', 'diaChXPred', '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=self.row_info_diamond['y_off'], srp=int(10000*self.row_info_diamond['pitch'])), temp_cuts)

    def Draw2DHistoDiamondCellOverlay(self, name, cells='all', cut=''):
        list_cuts = ['({l}<diaChYPred)&&(diaChYPred<{h})'.format(l=self.row_info_diamond['0'], h=self.row_info_diamond['0'] + self.row_info_diamond['pitch'] * self.row_info_diamond['num'])]
        if cells == 'good':
            list_cuts.append(self.goodAreasCutNames_diamond)
        elif cells == 'bad':
            list_cuts.append(self.badAreasCutNames_diamond)
        if cut != '':
            list_cuts.append(cut)
        temp_cuts = '&&'.join(list_cuts)
        self.Draw2DHisto(name, 0, self.col_pitch, self.cell_resolution, 'dia X [#mum]', 0, self.row_info_diamond['pitch'], self.cell_resolution, 'dia Y [#mum]', '((diaChXPred-{ox})*{p})%{p}'.format(ox=self.row_info_diamond['x_off'], p=self.col_pitch), '(((diaChYPred-{oy})*10000)%{srp})/10000'.format(oy=self.row_info_diamond['y_off'], srp=int(10000*self.row_info_diamond['pitch'])), temp_cuts)

    def DrawTCutCentersInCellOverlay(self, name):
        self.canvas[name].cd()
        self.tcutg_diamond_center.Draw('same')

    def FitLanGaus(self, name, conv_steps=100, color=ro.kRed):
        self.canvas[name].cd()
        self.langaus[name] = LanGaus(self.histo[name])
        self.langaus[name].LanGausFit(conv_steps)
        lowbin, highbin = self.histo[name].FindFirstBinAbove(0, 1), self.histo[name].FindLastBinAbove(0, 1)
        xlow, xhigh = self.histo[name].GetBinLowEdge(lowbin), self.histo[name].GetBinLowEdge(highbin + 1)
        self.line[name] = ro.TLine(xlow, 0, xhigh, 0)
        self.line[name].SetLineColor(ro.kViolet + 1)
        self.line[name].SetLineWidth(4)
        fitmean = self.langaus[name].fit.Mean(xlow, xhigh)
        self.langaus[name].fit.Draw('same')
        self.langaus[name].fit.SetLineColor(color)
        self.line[name].Draw('same')
        ro.gPad.Update()
        self.histo[name].FindObject('stats').SetOptFit(1)
        self.histo[name].FindObject('stats').SetX1NDC(0.6)
        self.histo[name].FindObject('stats').SetX2NDC(0.9)
        self.histo[name].FindObject('stats').SetY1NDC(0.6)
        self.histo[name].FindObject('stats').SetY2NDC(0.9)
        AddLineToStats(self.canvas[name], 'Mean_{Fit}', fitmean)
        self.histo[name].SetStats(0)
        self.canvas[name].Modified()
        ro.gPad.Update()
        print '{n}: <PH> = {f}'.format(n=name, f=fitmean)

    def Test(self, num=0, clust_size=2):
        self.cell_resolution = 50.0 / 13.0 if self.col_pitch == 50 else 100.0 / 51
        if num == 0:
            self.cell_resolution = 50.0 / 25.0
            self.ResetAreas()
            self.AddGoodAreasCol(2 + 1 - clust_size, 17, 20)
            self.AddGoodAreasCol(3 + 1 - clust_size, 7, 9)
            self.AddGoodAreasCol(3 + 1 - clust_size, 17, 19)
            self.AddGoodAreasCol(4 + 1 - clust_size, 6, 13)
            self.AddGoodAreasCol(4 + 1 - clust_size, 16, 20)
            # self.AddGoodAreasCol(5 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(5 + 1 - clust_size, 8, 16)
            self.AddGoodAreasCol(6 + 1 - clust_size, 7, 9)
            # self.AddGoodAreasCol(6 + 1 - clust_size, 11, 18)
            self.AddGoodAreasCol(6 + 1 - clust_size, 11, 17)
            # self.AddGoodAreasCol(7 + 1 - clust_size, 8, 8)
            # self.AddGoodAreasCol(7 + 1 - clust_size, 12, 12)
            self.AddGoodAreasCol(7 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(7 + 1 - clust_size, 21, 21)
            self.AddGoodAreasCol(8 + 1 - clust_size, 15, 21)
            self.AddGoodAreasCol(8 + 1 - clust_size, 24, 24)
            self.AddGoodAreasCol(9 + 1 - clust_size, 12, 12)
            self.AddGoodAreasCol(9 + 1 - clust_size, 14, 17)
            self.AddGoodAreasCol(9 + 1 - clust_size, 19, 25)
            self.AddGoodAreasCol(10 + 1 - clust_size, 12, 17)
            self.AddGoodAreasCol(10 + 1 - clust_size, 19, 19)
            self.AddGoodAreasCol(11 + 1 - clust_size, 13, 17)
            self.AddGoodAreasCol(11 + 1 - clust_size, 19, 20)
            self.RemoveFromGoodArea(3 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(6 + 1 - clust_size, 7)
            self.RemoveFromGoodArea(7 + 1 - clust_size, 14)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 16)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(10 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(10 + 1 - clust_size, 19)
            self.AddRemainingToBadAreas()
        elif num == 1:
            self.cell_resolution = 50.0 / 11.0
            # self.cell_resolution = 50.0 / 17.0
            self.ResetAreas()
            # self.AddGoodAreas(2 + 1 - clust_size, 8)
            self.AddGoodAreasCol(3 + 1 - clust_size, 7, 9)
            self.AddGoodAreasCol(4 + 1 - clust_size, 6, 13)
            self.AddGoodAreasCol(5 + 1 - clust_size, 8, 19)
            # self.AddGoodAreasCol(6 + 1 - clust_size, 11, 18)
            self.AddGoodAreasCol(6 + 1 - clust_size, 11, 16)
            # self.AddGoodAreasCol(7 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(7 + 1 - clust_size, 14, 17)
            # self.AddGoodAreasCol(8 + 1 - clust_size, 15, 21)
            self.AddGoodAreasCol(8 + 1 - clust_size, 15, 17)
            self.AddGoodAreasCol(9 + 1 - clust_size, 14, 17)
            self.AddGoodAreasCol(10 + 1 - clust_size, 12, 17)
            self.AddGoodAreasCol(11 + 1 - clust_size, 13, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 19)
            self.RemoveFromGoodArea(6 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(7 + 1 - clust_size, 14)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 16)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(10 + 1 - clust_size, 16)
            self.RemoveFromGoodArea(10 + 1 - clust_size, 17)
            self.AddRemainingToBadAreas()
        elif num == 2:
            self.cell_resolution = 50.0 / 15.0
            self.ResetAreas()
            self.AddGoodAreasCol(4 + 1 - clust_size, 6, 13)
            # self.AddGoodAreasCol(5 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(5 + 1 - clust_size, 8, 16)
            # self.AddGoodAreasCol(6 + 1 - clust_size, 11, 18)
            self.AddGoodAreasCol(6 + 1 - clust_size, 11, 17)
            # self.AddGoodAreasCol(7 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(7 + 1 - clust_size, 15, 19)
            self.AddGoodAreasCol(8 + 1 - clust_size, 15, 21)
            self.AddGoodAreasCol(8 + 1 - clust_size, 15, 21)
            self.AddRemainingToBadAreas()
        elif num == 3:
            self.cell_resolution = 50.0 / 15.0
            self.ResetAreas()
            self.AddGoodAreasCol(4 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(5 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(6 + 1 - clust_size, 8, 19)
            self.AddRemainingToBadAreas()
        elif num == 4:
            self.ResetAreas()
            # self.AddGoodAreasCol(5 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(5 + 1 - clust_size, 8, 16)
            # self.AddGoodAreasCol(6 + 1 - clust_size, 8, 19)
            self.AddGoodAreasCol(6 + 1 - clust_size, 8, 16)
            self.AddRemainingToBadAreas()
        elif num == 5:
            self.ResetAreas()
            self.AddGoodAreasCol(5 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(6 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(7 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(8 + 1 - clust_size, 14, 19)
            self.AddGoodAreasCol(9 + 1 - clust_size, 14, 19)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(6 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(7 + 1 - clust_size, 14)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 16)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(9 + 1 - clust_size, 18)
            self.AddRemainingToBadAreas()
        elif num == 6:
            self.cell_resolution = 50.0 / 11.0
            self.ResetAreas()
            self.AddGoodAreasCol(5 + 1 - clust_size, 14, 18)
            self.AddGoodAreasCol(6 + 1 - clust_size, 14, 18)
            self.AddGoodAreasCol(7 + 1 - clust_size, 14, 18)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 17)
            self.RemoveFromGoodArea(5 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(6 + 1 - clust_size, 18)
            self.RemoveFromGoodArea(7 + 1 - clust_size, 14)
            self.AddRemainingToBadAreas()
        elif num == 7:
            self.cell_resolution = 50.0 / 11.0
            self.ResetAreas()
            # self.AddGoodAreasCol(5 + 1 - clust_size, 15, 18)
            self.AddGoodAreasCol(5 + 1 - clust_size, 15, 16)
            # self.AddGoodAreasCol(6 + 1 - clust_size, 15, 18)
            self.AddGoodAreasCol(6 + 1 - clust_size, 15, 17)
            self.AddGoodAreasCol(7 + 1 - clust_size, 15, 18)
            self.AddGoodAreasCol(8 + 1 - clust_size, 15, 18)
            self.AddRemainingToBadAreas()
        elif num == 8:
            self.cell_resolution = 50.0 / 11.0
            self.ResetAreas()
            self.AddGoodAreasRow(15, 5 + 1 - clust_size, 11 + 2 - 2 * clust_size)
            self.AddGoodAreasRow(16, 5 + 1 - clust_size, 11 + 2 - 2 * clust_size)
            self.AddGoodAreasRow(17, 5 + 1 - clust_size, 11 + 2 - 2 * clust_size)
            self.AddRemainingToBadAreas()
        elif num == 9:
            self.cell_resolution = 2
            self.ResetAreas()
            self.SelectGoodAndBadByThreshold(800)
        elif num == 100:
            self.cell_resolution = 2
            self.ResetAreas()
            self.SelectGoodAndBadByThreshold(200)
        self.PlotTests(num, clust_size)

    def PlotTests(self, num=0, cluster_size=2):
        w = 0
        for ch in xrange(1, cluster_size + 1):
            #  plot map
            self.DrawProfile2DDiamond('ph{c}_map_test{n}'.format(c=ch, n=num), varz='clusterCharge'+str(ch), cuts='({l}<=diaChYPred)&&(diaChYPred<={h})'.format(l=tg.row_info_diamond['0'], h=tg.row_info_diamond['up']))
            # ro.gPad.Update()
            self.canvas['ph{c}_map_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.gridTextDiamond.Draw('same TEXT0')
            self.profile['ph{c}_map_test'.format(c=ch)+str(num)].GetXaxis().SetRangeUser(self.ch_ini - 1, self.ch_end + 1)
            self.profile['ph{c}_map_test'.format(c=ch)+str(num)].GetYaxis().SetRangeUser(self.row_info_diamond['0'] - int(self.row_info_diamond['pitch'] / self.bins_per_ch_y), self.row_info_diamond['up'] + int(self.row_info_diamond['pitch'] / self.bins_per_ch_y))
            #  draw selected areas
            self.DrawGoodAreasDiamond('ph{c}_map_test{n}'.format(c=ch, n=num))
            self.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}'.format(c=ch, n=num))
            #  plot cell overlay
            self.Draw2DProfileDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num), var='clusterCharge'+str(ch), cells='good')
            self.canvas['ph{c}_cells_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            #  show center region in cell
            self.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num))
            #  draw ph of selected areas
            self.DrawPHGoodAreas('ph{c}_test{n}'.format(c=ch, n=num), var='clusterCharge'+str(ch))
            self.canvas['ph{c}_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=ch, n=num), cells='good', var='clusterCharge'+str(ch))
            self.canvas['ph{c}_test{n}_centers'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            #  fit distribution for central region
            self.FitLanGaus('ph{c}_test{n}_centers'.format(c=ch, n=num), color=ro.kRed)
            #  get difference between cell and center
            self.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=ch, n=num))
            self.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Reset('ICES')
            self.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Add(self.histo['ph{c}_test{n}'.format(c=ch, n=num)], self.histo['ph{c}_test{n}_centers'.format(c=ch, n=num)], 1, -1)
            self.canvas['ph{c}_test{n}_periphery'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.FitLanGaus('ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlue)
            self.canvas['ph{c}_test{n}'.format(c=ch, n=num)].cd()
            self.langaus['ph{c}_test{n}_centers'.format(c=ch, n=num)].fit.Draw('same')
            self.langaus['ph{c}_test{n}_periphery'.format(c=ch, n=num)].fit.Draw('same')
            self.DrawDoubleLangaus('ph{c}_test{n}'.format(c=ch, n=num), 'ph{c}_test{n}_centers'.format(c=ch, n=num), 'ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlack)
            # ro.gPad.Update()
            #  position of negative clusters
            self.Draw2DHistoDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
            self.histo['hm{c}_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.ch_ini - 1, self.ch_end + 1)
            self.histo['hm{c}_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(self.row_info_diamond['0'] - int(self.row_info_diamond['pitch'] / self.bins_per_ch_y), self.row_info_diamond['up'] + int(self.row_info_diamond['pitch'] / self.bins_per_ch_y))
            self.canvas['hm{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawGoodAreasDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
            self.DrawBadAreasDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
            self.DrawGoodAreasDiamondCenters('hm{c}_test{n}_negative'.format(c=ch, n=num))
            self.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num), varz='clusterCharge'+str(ch), cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c}))'.format(c=self.neg_cut))
            self.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.ch_ini - 1, self.ch_end + 1)
            self.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(self.row_info_diamond['0'] - int(self.row_info_diamond['pitch'] / self.bins_per_ch_y), self.row_info_diamond['up'] + int(self.row_info_diamond['pitch'] / self.bins_per_ch_y))
            self.canvas['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawGoodAreasDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
            self.DrawBadAreasDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
            self.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
            self.Draw2DProfileDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=ch, n=num), var='clusterCharge'+str(ch), cells='good', cut='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=ch, n=num), var='clusterCharge'+str(ch), cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.FitLanGaus('ph{c}_test{n}_negative'.format(c=ch, n=num), color=ro.kBlue)
            self.Draw2DProfileDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge'+str(ch), cells='good', cut='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge'+str(ch), cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
            w += 15
            self.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=ch, n=num), color=ro.kRed)
        if cluster_size == 2:
            self.Draw2DProfileDiamondCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1', cells='good')
            self.canvas['ph2_minus_ph1_map_test{n}'.format(n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawTCutCentersInCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num))
            self.phmin, self.phmax = -1500, 1500
            self.phbins = int(np.floor((self.phmax - self.phmin) / 20.0 + 0.5))
            self.DrawPHGoodAreas('ph2_minus_ph1_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1')
            self.canvas['ph2_minus_ph1_test{n}'.format(n=num)].SetWindowPosition(w, w)
            w += 15
            # self.FitLanGaus('ph2_minus_ph1_test{n}'.format(n=num), color=ro.kViolet)
            self.phmin, self.phmax = 0, 4000
            self.phbins = 200
            self.Draw2DProfileDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=3,n=num), var='clusterChargeN', cells='good', cut='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_cells_test{n}_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=3, n=num), var='clusterChargeN', cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_test{n}_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
            self.FitLanGaus('ph{c}_test{n}_negative'.format(c=3, n=num), color=ro.kBlue)
            w += 15
            self.Draw2DProfileDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=3,n=num), var='clusterChargeN', cells='good', cut='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_cells_test{n}_no_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
            w += 15
            self.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=3, n=num), var='clusterChargeN', cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(c=self.neg_cut))
            self.canvas['ph{c}_test{n}_no_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
            w += 15
            self.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=3, n=num), color=ro.kRed)
            # ro.gPad.Update()
            # self.Draw2DProfileDiamondCellOverlay('ph_cells_not_test{n}'.format(n=num), cells='bad')
            # self.DrawPHBadAreas('ph_not_test{n}'.format(n=num))

    def SaveCanvasInlist(self, list, subdir):
        for canvas in list:
            if self.canvas.has_key(canvas):
                self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.png'.format(d=self.dir, r=self.run, sd=subdir, c=canvas))
                self.canvas[canvas].SaveAs('{d}/{r}/{sd}/{c}.root'.format(d=self.dir, r=self.run, sd=subdir, c=canvas))

    def DrawDoubleLangaus(self, name, name1, name2, color=ro.kBlack):
        if self.langaus.has_key(name1) and self.langaus.has_key(name2):
            langaus1 = tg.langaus[name1].fit
            langaus2 = tg.langaus[name2].fit
            self.doubleLangaus[name] = ro.TF1(name, self.TwoLanGaus, 0, 4000, 8)
            params1 = np.zeros(4, 'f8')
            params2 = np.zeros(4, 'f8')
            langaus1.GetParameters(params1)
            langaus2.GetParameters(params2)
            params12 = np.concatenate((params1, params2))
            self.doubleLangaus[name].SetNpx(1000)
            self.doubleLangaus[name].SetParameters(params12)
            self.doubleLangaus[name].SetParNames('Width1', 'MP1', 'Area1', 'GSigma1', 'Width2', 'MP2', 'Area2', 'GSigma2')
            fitmean = self.doubleLangaus[name].Mean(0, 4000)
            self.canvas[name].cd()
            self.doubleLangaus[name].Draw('same')
            self.doubleLangaus[name].SetLineColor(color)
            AddLineToStats(self.canvas[name], 'Mean_{Fit}', fitmean)
            self.histo[name].SetStats(0)
            self.canvas[name].Modified()
            ro.gPad.Update()
            print '{n}: <PH> = {f}'.format(n=name, f=fitmean)

    def TwoLanGaus(self, x, params):
        mpc1 = params[1] - self.mpshift * params[0]
        mpc2 = params[5] - self.mpshift * params[4]
        xlow1, xhigh1 = [x[0] + self.sigma_conv * i * params[3] for i in [-1, 1]]
        xlow2, xhigh2 = [x[0] + self.sigma_conv * i * params[7] for i in [-1, 1]]
        step1 = (xhigh1 - xlow1) / self.conv_steps
        step2 = (xhigh2 - xlow2) / self.conv_steps
        sums1 = 0
        sums2 = 0
        for i in xrange(1, int(np.ceil(self.conv_steps / 2.0 + 1))):
            xx1 = xlow1 + (i - 0.5) * step1
            xx2 = xlow2 + (i - 0.5) * step2
            fland1 = ro.TMath.Landau(xx1, mpc1, params[0]) / params[0]
            fland2 = ro.TMath.Landau(xx2, mpc2, params[4]) / params[4]
            sums1 += fland1 * ro.TMath.Gaus(x[0], xx1, params[3])
            sums2 += fland2 * ro.TMath.Gaus(x[0], xx2, params[7])
            xx1 = xhigh1 - (i - 0.5) * step1
            xx2 = xhigh2 - (i - 0.5) * step2
            fland1 = ro.TMath.Landau(xx1, mpc1, params[0]) / params[0]
            fland2 = ro.TMath.Landau(xx2, mpc2, params[4]) / params[4]
            sums1 += fland1 * ro.TMath.Gaus(x[0], xx1, params[3])
            sums2 += fland2 * ro.TMath.Gaus(x[0], xx2, params[7])
        return params[2] * step1 * sums1 / (np.sqrt(2 * np.pi, dtype='f8') * params[3]) + params[6] * step2 * sums2 / (np.sqrt(2 * np.pi, dtype='f8') * params[7])


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
    parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
    parser.add_option('-a', '--al', dest='al', action='store_true', default=True, help='enable find grid through data and alignment')
    parser.add_option('-c', '--cellsize', dest='cellsize', type='int', default=50, help='cell size of the square 3D device')
    parser.add_option('-t', '--test', dest='testnumber', type='int', default=-1, help='Run a automatically one of the predefined tests')

    (options, args) = parser.parse_args()
    run = int(options.run)
    dir = str(options.dir)
    use_align = bool(options.al)
    testnum = int(options.testnumber)
    cellsize = int(options.cellsize) if testnum != 100 else 100

    tg = TransparentGrid(dir=dir, run=run, cellsize=cellsize)
    tg.CheckFoldersAndFiles()
    tg.OpenFileAndGetTree()
    tg.FindDiamondChannelLimits()
    tg.SetLines(try_align=use_align)
    if testnum in (range(0, 10) + [100]):
        tg.CreateTCutGs()
        tg.Test(testnum)
        if not os.path.isdir('{d}/{r}/test{n}'.format(d=tg.dir, r=tg.run, n=testnum)):
            os.makedirs('{d}/{r}/test{n}'.format(d=tg.dir, r=tg.run, n=testnum))
        list_of_canvas = tg.canvas.keys()
        # list_of_canvas = [['ph{c}_map_test{n}'.format(c=ch, n=testnum), 'ph{c}_cells_test{n}'.format(c=ch, n=testnum), 'ph{c}_test{n}'.format(c=ch, n=testnum), 'ph{c}_test{n}_centers'.format(c=ch, n=testnum),
        #                   'ph{c}_test{n}_periphery'.format(c=ch, n=testnum), 'ph{c}_test{n}'.format(c=ch, n=testnum), 'ph2_minus_ph1_map_test{n}'.format(n=testnum), 'ph2_minus_ph1_test{n}'.format(n=testnum)] for ch in [1, 2]]
        # list_of_canvas = [elem for it in list_of_canvas for elem in it]
        tg.SaveCanvasInlist(list=list_of_canvas, subdir='test{n}'.format(n=testnum))
