#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
from copy import deepcopy
from NoiseExtraction import NoiseExtraction
import os

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64', 'ULong64_t': 'uint64', 'Bool_t': 'bool'}
diaChs = 128

class PedestalAnalysis:
    def __init__(self, run=22011, dir=''):
        print 'Creating NoiseExtraction instance for run:', run
        self.run = run
        self.dir = dir
        self.rootFile = self.pedTree = None
        self.adc_vect = self.ped_vect = self.sigma_vect = self.cm_vect = self.ped_cmc_vect = self.sigma_cmc_vect = None
        self.adc_hist, self.ped_hist, self.sigma_hist, self.cm_hist, self.ped_cmc_hist, self.sigma_cmc_hist = ro.TH2D(), ro.TH2D(), ro.TH2D(), ro.TH1D(), ro.TH2D(), ro.TH2D()
        self.dicBraNames = {'rawTree.DiaADC': 'adc_'+str(self.run), 'diaPedestalMean': 'ped_'+str(self.run), 'diaPedestalSigma': 'sigma_'+str(self.run), 'diaPedestalMeanCMN': 'ped_cmc_'+str(self.run), 'diaPedestaSigmaCMN': 'sig_cmc_'+str(self.run), 'commonModeNoise': 'cm_'+str(self.run)}
        self.dicBraVectChs = {}
        self.dicBraVect1ch = {}
        self.allBranches = self.dicBraNames.keys()
        for key in self.allBranches:
            if key != 'commonModeNoise':
                self.dicBraVectChs[key] = None
            else:
                self.dicBraVect1ch[key] = None
        self.entries = 0
        self.dictBraHist = {'rawTree.DiaADC': self.adc_hist, 'diaPedestalMean': self.ped_hist, 'diaPedestalSigma': self.sigma_hist, 'diaPedestalMeanCMN': self.ped_cmc_hist, 'diaPedestaSigmaCMN': self.sigma_cmc_hist}
        self.ev_axis = {'min':0, 'max': 0, 'bins': 1000}
        self.ch_axis = {'min': -0.5, 'max': diaChs - 0.5, 'bins': diaChs}
        # self.channels = range(28, 48) if connect==1 else range(75, 113) if connect == 0 else range(chL, chH + 1)
        # self.sourceDir = sourceDir
        # self.outputDir = outputDir
        # self.noise = {ch: NoiseExtraction(self.run, self.sourceDir, self.outputDir) for ch in self.channels}
        # self.noiseAll = NoiseExtraction(self.run, self.sourceDir, self.outputDir)
        # for ch in self.channels:
        #     self.noise[ch].SetChannelsBatch(str(ch), ch, ch)
        #     self.noise[ch].Get2DMap(True)
        #     self.noise[ch].Get1DHistos(True)
        # self.noiseAll.SetChannelsBatch('{l}_{h}'.format(l=min(self.channels), h=max(self.channels)), int(min(self.channels)), int(max(self.channels)))
        # self.noiseAll.Get2DMap(True)
        # self.noiseAll.Get1DHistos(True)
        # self.fitNoiseCMNC = None
        # self.fitNoise = None
        # self.bla = []
        # self.noiseMeans = [self.noise[ch].noiseHisto1D.GetMean() for ch in self.channels]
        # self.noiseCMNCMeans = [self.noise[ch].noiseCMNHisto1D.GetMean() for ch in self.channels]
        # self.noiseSigma = [self.noise[ch].noiseHisto1D.GetRMS() for ch in self.channels]
        # self.noiseCMNCSigma = [self.noise[ch].noiseCMNHisto1D.GetRMS() for ch in self.channels]
        # print 'The mean of the channels\' mean before CMNC is: {m} +/- {s}'.format(m=average(self.noiseMeans), s=std(self.noiseMeans))
        # print 'The mean of the channels\' sigma before CMNC is: {m} +/- {s}'.format(m=average(self.noiseSigma), s=std(self.noiseSigma))
        # print 'The mean of the channels\' mean after CMNC is: {m} +/- {s}'.format(m=average(self.noiseCMNCMeans), s=std(self.noiseCMNCMeans))
        # print 'The mean of the channels\' sigma after CMNC is: {m} +/- {s}'.format(m=average(self.noiseCMNCSigma), s=std(self.noiseCMNCSigma))
        # print 'The mean of all the channels\' mean before CMNC is: {m} +/- {s}'.format(m=self.noiseAll.noiseHisto1D.GetMean(), s=self.noiseAll.noiseHisto1D.GetRMS())
        # print 'The mean of all the channels\' mean after CMNC is: {m} +/- {s}'.format(m=self.noiseAll.noiseCMNHisto1D.GetMean(), s=self.noiseAll.noiseCMNHisto1D.GetRMS())

    def LoadROOTFile(self):
        self.rootFile = ro.TFile('{d}/pedestalData.{r}.root'.format(d=self.dir, r=self.run), 'READ')
        self.pedTree = self.rootFile.Get('pedestalTree')
        self.entries = self.pedTree.GetEntries()

    def LoadVectorsFromBranches(self, first_ev=0):
        if self.pedTree is None:
            self.LoadROOTFile()

        num_bra_chs = len(self.dicBraVectChs)
        if num_bra_chs < 1:
            print 'the dictionary of branches and vectors is empty! try again'
            return
        channels = self.pedTree.GetLeaf(self.dicBraVectChs[self.dicBraVectChs.keys()[0]]).GetLen()
        for branch in self.dicBraVectChs.iterkeys():
            if self.pedTree.GetLeaf(branch).GetLen() != channels:
                print 'Given branches have different sizes! try again'
                return
        branchesChs = self.dicBraVectChs.keys()
        leng = self.pedTree.Draw(':'.join(branchesChs), '', 'goff para', self.entries, first_ev)
        if leng == -1:
            print 'Error, could not load the branches. try again'
            return
        while leng > self.pedTree.GetEstimate():
            self.pedTree.SetEstimate(leng)
            leng = self.pedTree.Draw(':'.join(branchesChs), '', 'goff para', self.entries, first_ev)
        self.entries = leng / channels
        for pos, branch in enumerate(branchesChs):
            type = dicTypes[self.pedTree.GetLeaf(branch).GetTypeName()]
            self.dicBraVectChs[branch] = self.pedTree.GetVal(pos)
            self.dicBraVectChs[branch] = np.array([[self.dicBraVectChs[branch][ev * channels + ch] for ch in xrange(channels)] for ev in xrange(self.entries)], dtype=np.dtype(type))

        num_bra_1ch = len(self.dicBraVect1ch)
        if num_bra_1ch < 1:
            print 'the dictionary of branches and vectors is empty! try again'
            return
        channel = 1
        for branch in self.dicBraVect1ch.iterkeys():
            if self.pedTree.GetLeaf(branch).GetLen() != channels:
                print 'Given branches have different sizes different to 1! try again'
                return
        branches1ch = self.dicBraVect1ch.keys()
        leng = self.pedTree.Draw(':'.join(branches1ch), '', 'goff para', self.entries, first_ev)
        if leng == -1:
            print 'Error, could not load the branches. try again'
            return
        while leng > self.pedTree.GetEstimate():
            self.pedTree.SetEstimate(leng)
            leng = self.pedTree.Draw(':'.join(branches1ch), '', 'goff para', self.entries, first_ev)
        for pos, branch in enumerate(branches1ch):
            type = dicTypes[self.pedTree.GetLeaf(branch).GetTypeName()]
            self.dicBraVect1ch[branch] = self.pedTree.GetVal(pos)
            self.dicBraVect1ch[branch] = np.array([self.dicBraVect1ch[branch][ev] for ev in xrange(self.entries)], dtype=np.dtype(type))

        self.adc_vect, self.ped_vect, self.sigma_vect, self.ped_cmc_vect, self.sigma_cmc_vect = self.dicBraVectChs['rawTree.DiaADC'], self.dicBraVectChs['diaPedestalMean'], self.dicBraVectChs['diaPedestalSigma'], self.dicBraVectChs['diaPedestalMeanCMN'], self.dicBraVectChs['diaPedestaSigmaCMN']
        self.cm_vect = self.dicBraVect1ch['commonModeNoise']


    def CreateHistogramsForVectors(self):
        for branch in self.allBranches:
            name = self.dicBraNames[branch]
            if branch in self.dicBraVect1ch.keys():
                ymin, ymax = self.dicBraVect1ch[branch].min(), self.dicBraVect1ch[branch].max()
                ybins = 500.
                self.dictBraHist[branch].SetNameTitle(name, name)
                self.dictBraHist[branch].GetXaxis().Set(self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'])
                self.dictBraHist[branch].GetYaxis().Set(ybins, ymin, ymax)
            elif branch in self.dicBraVectChs.keys():
                zmin, zmax = self.dicBraVectChs[branch].min(), self.dicBraVectChs[branch].max()
                zbins = 500.
                self.dictBraHist[branch].SetNameTitle(name, name)
                self.dictBraHist[branch].GetXaxis().Set(self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'])
                self.dictBraHist[branch].GetYaxis().Set(self.ch_axis['bins'], self.ch_axis['min'], self.ch_axis['max'])
                self.dictBraHist[branch].GetZaxis().Set(zbins, zmin, zmax)

    def FillHistograms(self):
        for ev in xrange(self.entries):
            self.cm_hist.Fill(ev, self.cm_vect[ev])
            for ch in xrange(diaChs):
                self.adc_hist.Fill(ev, ch, self.adc_vect[ev, ch])
                self.ped_hist.Fill(ev, ch, self.ped_vect[ev, ch])
                self.sigma_hist.Fill(ev, ch, self.sigma_vect[ev, ch])
                self.ped_cmc_hist.Fill(ev, ch, self.ped_cmc_vect[ev, ch])
                self.sigma_cmc_hist.Fill(ev, ch, self.sigma_cmc_vect[ev, ch])

                
                # def Get1DHistos(self):
                    #     if self.fidcut == 0 or self.fidcutCMN == 0:
                    #         self.CreateFidCut()
                    #     self.noiseHisto1D = self.noiseHisto2D.ProjectionY('{n}_all_Chs_Region'.format(n=self.noiseHisto2D.GetName()), int(self.noiseHisto2D.GetXaxis().GetBinLowEdge(self.diamondNoiseChannels['ch_low'])), int(self.noiseHisto2D.GetXaxis().GetBinUpEdge(self.diamondNoiseChannels['ch_high'])+1))
                    #     self.noiseHisto1D.GetXaxis().SetRangeUser(-32, 32)
                    #     self.noiseHisto1D.SetFillColor(38)
                    #     canvas_name = 'c_noise_1D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
                    #     canvas = self.CreateCanvas(canvas_name)
                    #     gStyle.SetOptStat('n')
                    #     canvas.cd()
                    #     self.noiseHisto1D.SetLineWidth(3)
                    #     self.noiseHisto1D.GetYaxis().SetTitle('entries')
                    #     self.noiseHisto1D.GetYaxis().SetTitleOffset(1.5)
                    #     self.FitHistogramGaus(self.noiseHisto1D, self.fitNoise, '{n}_fit'.format(n=self.noiseHisto1D.GetName()))
                    #     gStyle.SetOptFit(0111)
                    #     self.noiseHisto1D.Draw()
                    #     name = 'histo1D_noise_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
                    #     self.SaveCanvas(canvas, name)
                    #     self.bla.append(canvas)
                    #
                    #     self.noiseCMNHisto1D = self.noiseCMNHisto2D.ProjectionY('{n}_all_Chs_Region'.format(n=self.noiseCMNHisto2D.GetName()),
                    #                                                       int(self.noiseCMNHisto2D.GetXaxis().GetBinLowEdge(
                    #                                                           self.diamondNoiseChannels['ch_low'])),
                    #                                                       int(self.noiseCMNHisto2D.GetXaxis().GetBinUpEdge(
                    #                                                           self.diamondNoiseChannels['ch_high'])+1))
                    #     self.noiseCMNHisto1D.GetXaxis().SetRangeUser(-32, 32)
                    #     self.noiseCMNHisto1D.SetFillColor(38)
                    #     canvas_nameCMN = 'c_noiseCMNC_1D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                    #                                                            l=self.diamondNoiseChannels['ch_low'],
                    #                                                            h=self.diamondNoiseChannels['ch_high'])
                    #     canvasCMN = self.CreateCanvas(canvas_nameCMN)
                    #     gStyle.SetOptStat('n')
                    #     canvasCMN.cd()
                    #     self.noiseCMNHisto1D.SetLineWidth(3)
                    #     self.noiseCMNHisto1D.GetYaxis().SetTitle('entries')
                    #     self.noiseCMNHisto1D.GetYaxis().SetTitleOffset(1.5)
                    #     self.FitHistogramGaus(self.noiseCMNHisto1D, self.fitNoiseCMNC, '{n}_fit'.format(n=self.noiseCMNHisto1D.GetName()))
                    #     gStyle.SetOptFit(0111)
                    #     self.noiseCMNHisto1D.Draw()
                    #     nameCMN = 'histo1D_noiseCMNC_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=self.diamondNoiseChannels['ch_low'], h=self.diamondNoiseChannels['ch_high'])
                    #     self.SaveCanvas(canvasCMN, nameCMN)
                    #     self.bla.append(canvasCMN)
                    #
                    # def FitHistogramGaus(self, histo, fit, fit_name):
                    #     mean = histo.GetMean()
                    #     rms = histo.GetRMS()
                    #     low = mean - 2*rms
                    #     up = mean + 2*rms
                    #     fit = TF1(fit_name, 'gaus', low, up)
                    #     histo.Fit(fit_name, 'rq')
                    #     fit.SetLineColor(kRed)
                    #     fit.SetLineWidth(5)
                    #
                    # def Get2DMap(self):
                    #     if self.fidcut == 0 or self.fidcutCMN == 0:
                    #         self.CreateFidCut()
                    #     self.noiseHisto2D.GetXaxis().SetRangeUser(0, diaChs)
                    #     self.noiseHisto2D.GetYaxis().SetRangeUser(-32, 32)
                    #     self.noiseHisto2D.GetZaxis().SetTitleOffset(1.2)
                    #     self.noiseCMNHisto2D.GetXaxis().SetRangeUser(0, diaChs)
                    #     self.noiseCMNHisto2D.GetYaxis().SetRangeUser(-32, 32)
                    #     self.noiseCMNHisto2D.GetZaxis().SetTitleOffset(1.2)
                    #
                    #     canvas_name = 'c_noise_2D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvas = self.CreateCanvas(canvas_name)
                    #     gStyle.SetOptStat('ne')
                    #     canvas.cd()
                    #     self.noiseHisto2D.Draw('colz')
                    #     self.fidcut.Draw('same')
                    #     self.bla.append(canvas)
                    #     name = 'histo2D_noise_nonhit_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvas, name)
                    #
                    #     canvasCMN_name = 'c_noiseCMNC_2D_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvasCMN = self.CreateCanvas(canvasCMN_name)
                    #     gStyle.SetOptStat('ne')
                    #     canvasCMN.cd()
                    #     self.noiseCMNHisto2D.Draw('colz')
                    #     self.fidcutCMN.Draw('same')
                    #     self.bla.append(canvasCMN)
                    #     nameCMN = 'histo2D_noiseCMNC_nonhit_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvasCMN, nameCMN)
                    #
                    # def Get2DMapFiducial(self):
                    #     if self.fidcut == 0 or self.fidcutCMN == 0:
                    #         self.CreateFidCut()
                    #     self.noiseHisto2D.GetXaxis().SetRangeUser(0, diaChs)
                    #     self.noiseHisto2D.GetYaxis().SetRangeUser(-32, 32)
                    #     self.noiseCMNHisto2D.GetXaxis().SetRangeUser(0, diaChs)
                    #     self.noiseCMNHisto2D.GetYaxis().SetRangeUser(-32, 32)
                    #
                    #     canvas_name = 'c_noise_2D_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvas = self.CreateCanvas(canvas_name)
                    #     gStyle.SetOptStat('ne')
                    #     canvas.cd()
                    #     self.noiseHisto2D.Draw('colz')
                    #     self.fidcut.Draw('same')
                    #     self.bla.append(canvas)
                    #     name = 'histo2D_noise_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvas, name)
                    #
                    #     canvasCMN_name = 'c_noiseCMNC_2D_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvasCMN = self.CreateCanvas(canvasCMN_name)
                    #     gStyle.SetOptStat('ne')
                    #     canvasCMN.cd()
                    #     self.noiseCMNHisto2D.Draw('colz')
                    #     self.fidcutCMN.Draw('same')
                    #     self.bla.append(canvasCMN)
                    #     nameCMN = 'histo2D_noiseCMNC_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection, l=int(self.diamondNoiseChannels['ch_low']), h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvasCMN, nameCMN)
                    #
                    #     self.noiseHisto2D_fid = self.noiseHisto2D.Clone('{n}_region'.format(n=self.noiseHisto2D.GetName()))
                    #     self.noiseHisto2D_fid.SetTitle(self.noiseHisto2D_fid.GetName())
                    #     nbinsMap = (self.noiseHisto2D.GetNbinsY() * self.noiseHisto2D.GetNbinsX() + 1)
                    #     for bin in xrange(1, nbinsMap):
                    #         x, y, z = Long(0), Long(0), Long(0)
                    #         if self.noiseHisto2D.GetBinContent(bin) > 0:
                    #             self.noiseHisto2D.GetBinXYZ(bin, x, y, z)
                    #             point = {'x': self.noiseHisto2D.GetXaxis().GetBinCenter(x), 'y': self.noiseHisto2D.GetYaxis().GetBinCenter(y)}
                    #             # if not self.WindingIsPointInPoly(point):
                    #             if not bool(int(self.fidcut.IsInside(point['x'], point['y']))):
                    #                 self.noiseHisto2D_fid.SetBinContent(bin, 0)
                    #     canvas_name = 'c_noise_2D_region_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                    #                                                          l=int(self.diamondNoiseChannels['ch_low']),
                    #                                                          h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvas = self.CreateCanvas(canvas_name)
                    #     gStyle.SetOptStat('n')
                    #     canvas.cd()
                    #     self.noiseHisto2D_fid.Draw('colz')
                    #     self.fidcut.Draw('same')
                    #     self.bla.append(canvas)
                    #     name = 'histo2D_noise_region_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                    #                                                             l=int(self.diamondNoiseChannels['ch_low']),
                    #                                                             h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvas, name)
                    #
                    #     self.noiseCMNHisto2D_fid = self.noiseCMNHisto2D.Clone('{n}_region'.format(n=self.noiseCMNHisto2D.GetName()))
                    #     self.noiseCMNHisto2D_fid.SetTitle(self.noiseCMNHisto2D_fid.GetName())
                    #     nbinsMap = (self.noiseCMNHisto2D.GetNbinsY() * self.noiseCMNHisto2D.GetNbinsX() + 1)
                    #     for bin in xrange(1, nbinsMap):
                    #         x, y, z = Long(0), Long(0), Long(0)
                    #         if self.noiseCMNHisto2D.GetBinContent(bin) > 0:
                    #             self.noiseCMNHisto2D.GetBinXYZ(bin, x, y, z)
                    #             point = {'x': self.noiseCMNHisto2D.GetXaxis().GetBinCenter(x),
                    #                      'y': self.noiseCMNHisto2D.GetYaxis().GetBinCenter(y)}
                    #             # if not self.WindingIsPointInPoly(point):
                    #             if not bool(int(self.fidcutCMN.IsInside(point['x'], point['y']))):
                    #                 self.noiseCMNHisto2D_fid.SetBinContent(bin, 0)
                    #     canvas_nameCMN = 'c_noiseCMNC_2D_region_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                    #                                                                 l=int(self.diamondNoiseChannels['ch_low']),
                    #                                                                 h=int(self.diamondNoiseChannels['ch_high']))
                    #     canvasCMN = self.CreateCanvas(canvas_nameCMN)
                    #     gStyle.SetOptStat('n')
                    #     canvasCMN.cd()
                    #     self.noiseCMNHisto2D_fid.Draw('colz')
                    #     self.fidcutCMN.Draw('same')
                    #     self.bla.append(canvasCMN)
                    #     nameCMN = 'histo2D_noiseCMNC_region_nonhit_fid_{r}_{n}_ch_[{l}_{h}]'.format(r=self.run, n=self.nameChSelection,
                    #                                                                    l=int(self.diamondNoiseChannels['ch_low']),
                    #                                                                    h=int(self.diamondNoiseChannels['ch_high']))
                    #     self.SaveCanvas(canvasCMN, nameCMN)
                    #
                    # def CreateFidCut(self):
                    #     if self.fidcut == 0 and self.fidcutCMN == 0:
                    #         self.fidcut = TCutG('fidcut', 5)
                    #         self.fidcut.SetVarX(self.noiseHisto2D.GetXaxis().GetName())
                    #         self.fidcut.SetVarY(self.noiseHisto2D.GetYaxis().GetName())
                    #         self.fidcut.SetPoint(0, self.diamondNoiseChannels['ch_low'], -32)
                    #         self.fidcut.SetPoint(1, self.diamondNoiseChannels['ch_low'], 32)
                    #         self.fidcut.SetPoint(2, self.diamondNoiseChannels['ch_high']+1, 32)
                    #         self.fidcut.SetPoint(3, self.diamondNoiseChannels['ch_high']+1, -32)
                    #         self.fidcut.SetPoint(4, self.diamondNoiseChannels['ch_low'], -32)
                    #         self.fidcutCMN = TCutG('fidcutCMN', 5)
                    #         self.fidcutCMN.SetVarX(self.noiseCMNHisto2D.GetXaxis().GetName())
                    #         self.fidcutCMN.SetVarY(self.noiseCMNHisto2D.GetYaxis().GetName())
                    #         self.fidcutCMN.SetPoint(0, self.diamondNoiseChannels['ch_low'], -32)
                    #         self.fidcutCMN.SetPoint(1, self.diamondNoiseChannels['ch_low'], 32)
                    #         self.fidcutCMN.SetPoint(2, self.diamondNoiseChannels['ch_high']+1, 32)
                    #         self.fidcutCMN.SetPoint(3, self.diamondNoiseChannels['ch_high']+1, -32)
                    #         self.fidcutCMN.SetPoint(4, self.diamondNoiseChannels['ch_low'], -32)
                    #     if self.fidcut != 0 and self.fidcutCMN != 0:
                    #         self.fidcut.SetLineColor(kRed)
                    #         self.fidcut.SetLineWidth(3)
                    #         self.fidcutCMN.SetLineColor(kRed)
                    #         self.fidcutCMN.SetLineWidth(3)
                    #
                    # def SetDefaultChannels(self):
                    #     self.diamondNoiseChannels = {'ch_low': 74, 'ch_high': 80}
                    #     self.nameChSelection = 'Strip'
                    #
                    # def SetChannels(self):
                    #     self.nameChSelection = raw_input('Enter the name for the channel selection: ')
                    #     low, high = '', ''
                    #     while True:
                    #         if low == '':
                    #             low = raw_input('Enter the lowest channel to analyse (int): ')
                    #         if not self.CheckStringInt(low):
                    #             print 'The channel entered is not an integer. Try again...'
                    #             low = raw_input('Enter the lowest channel to analyse (int): ')
                    #         elif 0 <= int(low) <= 127:
                    #             low = int(low)
                    #             break
                    #         else:
                    #             print 'The channel must be between [0, 127]. Try again...'
                    #             low = raw_input('Enter the lowest channel to analyse (int): ')
                    #     while True:
                    #         if high == '':
                    #             high = raw_input('Enter the highest channel to analyse (int): ')
                    #         if not self.CheckStringInt(high):
                    #             print 'The channel entered is not an integer. Try again...'
                    #             high = raw_input('Enter the highest channel to analyse (int): ')
                    #         elif low <= int(high) <= 127:
                    #             high = int(high)
                    #             break
                    #         else:
                    #             print 'The channel must be between [{l}, 127]. Try again...'.format(l=low)
                    #             high = raw_input('Enter the highest channel to analyse (int): ')
                    #     self.diamondNoiseChannels = {'ch_low': low, 'ch_high': high}
                    #     self.CreateFidCut()
                    #
                    # def CheckStringInt(self, s):
                    #     if s[0] in ('-', '+'):
                    #         return s[1:].isdigit()
                    #     return s.isdigit()
                    #
                    # def CreateCanvas(self, name='c0', w=1024, h=1024):
                    #     c = TCanvas(name, name, w, h)
                    #     c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
                    #     return deepcopy(c)
                    #
                    # def SaveCanvas(self, canvas, name):
                    #     canvas.SaveAs(self.outputDir + '/root/{n}.root'.format(n=name))
                    #     canvas.SaveAs(self.outputDir + '/Plots/{n}.png'.format(n=name))

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
    parser.add_option('-d', '--dir', dest='dir', default='.', type='string', help='source folder containing processed data of different runs')
    # parser.add_option('-o', '--outputDir', dest='output', default='', type='string', help='output folder containing the analysed results')
    # parser.add_option('-c', '--connected', dest='connect', default=1, type='int', help='do connected channels (1) or disconnected (0). Other integer value you would have to specify the range using -l and -H')
    # parser.add_option('-l', '--lowChannel', dest='low', default=0, type='int', help='lower channel to analyse e.g. 1')
    # parser.add_option('-H', '--highChannel', dest='high', default=0, type='int', help='higher channel to analyse e.g. 2')

    (options, args) = parser.parse_args()
    run = int(options.run)
    dir = str(options.dir)
    # output = str(options.output)
    # connect = int(options.connect)
    # low = int(options.low)
    # high = int(options.high)

    z = PedestalAnalysis(run, dir)
