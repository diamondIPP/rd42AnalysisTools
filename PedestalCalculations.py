# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
# sys.path.append('/home/sandiego/rd42AnalysisTools')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
# from Settings import Settings
from Utils import *
from PedestalDeviceCalculations import PedestalDeviceCalculations
import time
import numpy as np
import subprocess as subp
# from threading import Thread
# from multiprocessing import Manager, sharedctypes
import multiprocessing as mp
import multiprocessing.sharedctypes as mpsc
import ctypes
import ipdb

__author__ = 'DA'

class PedestalCalculations(mp.Process):
    def __init__(self, origTree, outdir, run_no, sliding_length=50, hit_fact=3, seed_fact=5, ev_ini=0, ev_end=0, showprogress=True):
        mp.Process.__init__(self)
        print 'Creating PedestalCalculations instance'
        self.rootTree = origTree
        self.out_dir = outdir
        self.run_no = run_no
        self.events = self.rootTree.GetEntries()
        self.raw_dut_branch = 'DiaADC'
        self.chs = self.rootTree.GetLeaf(self.raw_dut_branch).GetLen()
        self.dut_np_data_type, self.dut_root_data_type = 'uint16', 's'
        self.hit_fact = hit_fact
        self.seed_fact = seed_fact
        self.ev_ini = int(ev_ini)
        self.ev_end = int(ev_end)
        self.showprogress = showprogress

        # self.out_dir = self.settings.output_dir
        # self.sub_dir = self.settings.sub_dir
        # self.file_name = self.settings.file_name
        # self.run = self.settings.run
        # self.tree_name = self.settings.tree_name
        self.num_parallel = RoundInt(mp.cpu_count() / 2.0)
        self.slide_leng = sliding_length
        self.cm_cut = 4
        # self.ana_dir = self.settings.analysis_path
        self.ped_branches = ['diaChSignal', 'diaChSignalNoCmc', 'diaChHits', 'diaChSeed', 'diaChPed', 'diaChPedMeanCmc', 'diaChPedSigmaCmc', 'diaChCm', 'cmn', 'diaChHitsNoCmc', 'diaChSeedNoCmc', 'diaChPedNoCmc', 'diaChPedMeanNoCmc', 'diaChPedSigmaNoCmc']
        self.devices = ['dut']
        self.device_ped_process = {}
        self.outFile = ro.TFile()
        self.outTree = ro.TTree()

        self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs = np.zeros(self.chs, '?'), np.zeros(self.chs, '?'), np.zeros(self.chs, '?')
        self.GetScreenedChannels()

        # self.dic_device_to_pos = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}

        # Temporary numpy arrays to create the ctype vectors for multiprocessing with shared memory
        dut_ADC_signal_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_mean_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_sigma_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_is_ped_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_hit_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_seed_all = np.zeros((self.chs, self.events), dtype='uint8')

        # Temporary numpy arrays to create the ctype vectors for multiprocessing with shared memory for cmc calculations
        dut_ADC_signal_cmc_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_chs_cm_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_cm_all = np.zeros(self.events, dtype='float32')
        dut_ADC_mean_cmc_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_sigma_cmc_all = np.zeros((self.chs, self.events), dtype='float32')
        dut_ADC_is_ped_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_hit_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')
        dut_ADC_is_seed_cmc_all = np.zeros((self.chs, self.events), dtype='uint8')

        # The ctype vectors that will be shared in the multiprocessing
        # self.dut_ADC_signal_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_signal_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_signal_all))
        self.dut_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_mean_all))
        self.dut_ADC_sigma_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_sigma_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_sigma_all))
        self.dut_ADC_is_ped_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_ped_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_ped_all))
        self.dut_ADC_is_hit_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_hit_all))
        self.dut_ADC_is_seed_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_seed_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_seed_all))
        del dut_ADC_signal_all, dut_ADC_mean_all, dut_ADC_sigma_all, dut_ADC_is_ped_all, dut_ADC_is_hit_all, dut_ADC_is_seed_all

        # The ctype vectors that will be shared in the multiprocessing for cmc calculations
        # self.dut_ADC_signal_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_signal_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_signal_cmc_all))
        self.dut_ADC_chs_cm_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_chs_cm_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_chs_cm_all))
        self.dut_ADC_cm_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_cm_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_cm_all))
        self.dut_ADC_mean_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_mean_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_mean_cmc_all))
        self.dut_ADC_sigma_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_sigma_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_sigma_cmc_all))
        self.dut_ADC_is_ped_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_ped_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_ped_cmc_all))
        self.dut_ADC_is_hit_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_hit_cmc_all))
        self.dut_ADC_is_seed_cmc_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_seed_cmc_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_seed_cmc_all))
        del dut_ADC_signal_cmc_all, dut_ADC_chs_cm_all, dut_ADC_cm_all, dut_ADC_mean_cmc_all, dut_ADC_sigma_cmc_all, dut_ADC_is_ped_cmc_all, dut_ADC_is_hit_cmc_all, dut_ADC_is_seed_cmc_all

        # ctype vectors with all the ADC events for telescope and dut. The each process will read from this vectors
        self.dut_ADC_all_mp = self.LoadDutADCs()

        # numpy arrays to fill tree
        self.dut_ADC_all = np.zeros(self.chs, self.dut_np_data_type)
        self.dut_ADC_is_ped = np.zeros(self.chs, '?')
        self.dut_ADC_is_hit = np.zeros(self.chs, '?')
        self.dut_ADC_is_seed = np.zeros(self.chs, '?')
        self.dut_ADC_signal = np.zeros(self.chs, 'float32')
        self.dut_ADC_mean = np.zeros(self.chs, 'float32')
        self.dut_ADC_sigma = np.zeros(self.chs, 'float32')
        self.dut_ADC_cm_chs = np.zeros(self.chs, '?')
        self.dut_ADC_cm = np.zeros(1, 'float32')
        self.dut_ADC_is_ped_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_is_hit_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_is_seed_cmc = np.zeros(self.chs, '?')
        self.dut_ADC_signal_cmc = np.zeros(self.chs, 'float32')
        self.dut_ADC_mean_cmc = np.zeros(self.chs, 'float32')
        self.dut_ADC_sigma_cmc = np.zeros(self.chs, 'float32')
        self.fill_slide_leng = np.zeros(1, 'uint16')
        self.fill_slide_leng.fill(self.slide_leng)

        # branches to fill tree
        # self.b1, self.b2, self.b3, self.b4, self.b5, self.b6, self.b7, self.b8, self.b9, self.b10, self.b11, self.b12, self.b13, self.b14, self.b15, self.b16, self.b17 = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

    def GetScreenedChannels(self):
        list_bra = ['diaNcChs', 'diaNoisyChs', 'diaMaskedChs']
        Draw_Branches_For_Get_Val(self.rootTree, list_bra)
        Get_Branches_Value_To_Numpy(self.rootTree, list_bra, [self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs], 1, self.chs)

    def LoadDutADCs(self):
        # self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        # self.rootTree.SetEstimate(256 * self.rootTree.GetEntries())
        time0 = time.time()
        if self.showprogress:
            print 'Getting all events for DUT...', ; sys.stdout.flush()

        Draw_Branches_For_Get_Val(self.rootTree, [self.raw_dut_branch], 0, self.events, option='goff')
        temp = np.zeros((self.chs, self.events), dtype=self.dut_np_data_type)
        Get_Branches_Value_To_Numpy(self.rootTree, [self.raw_dut_branch], [temp], self.events, self.chs)
        if self.showprogress: print 'Done in', time.time() - time0, 'seconds'
        return mp.Array(np.ctypeslib.as_ctypes(temp)._type_, np.ctypeslib.as_ctypes(temp))

    def run(self):
        self.CalculateDevicesPedestals()

    def CalculateDevicesPedestals(self):
        jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
        for row in jobs_chunks:
            self.device_ped_process = []
            for it, device in enumerate(row):
                if (it == len(row) - 1) and self.showprogress:
                    print 'Calculating pedestals for', device, '. The progress is shown below:'
                    temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, self.dut_ADC_signal_all_mp, self.dut_ADC_signal_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)

                else:
                    # print 'Calculating pedestals for', device
                    temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations(self.out_dir, self.run_no, self.hit_fact, self.seed_fact, self.dut_np_data_type, self.chs, self.events, self.diaNcChs, self.diaNoisyChs, self.diaMaskedChs, self.slide_leng, self.cm_cut, device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, self.dut_ADC_signal_all_mp, self.dut_ADC_signal_cmc_all_mp, -1)
                    # temp = PedestalDeviceCalculations('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_ped_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, self.dut_ADC_chs_cm_all_mp, self.dut_ADC_cm_all_mp, self.dut_ADC_mean_cmc_all_mp, self.dut_ADC_sigma_cmc_all_mp, self.dut_ADC_is_ped_cmc_all_mp, self.dut_ADC_is_hit_cmc_all_mp, self.dut_ADC_is_seed_cmc_all_mp, -1)
                temp.start()
                self.device_ped_process.append(temp)
            for j in self.device_ped_process:
                j.join()

        # self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='UPDATE')
        self.outFile = ro.TFile('{d}/{r}/pedestal.{s}.{r}.root'.format(d=self.out_dir, r=self.run_no, s=self.slide_leng), 'RECREATE')
        self.outTree = ro.TTree('pedTree', 'pedTree')
        self.SetBranches()
        self.FillTree()


    def SetBranches(self):
        if self.showprogress:
            print 'Setting branches in tree...', ; sys.stdout.flush()
        # self.outTree.Branch('diaChPedNoCmc', self.dut_ADC_is_ped, 'diaChPedNoCmc[{chs}]/O'.format(chs=self.chs))
        # self.outTree.Branch('diaChHitsNoCmc', self.dut_ADC_is_hit, 'diaChHitsNoCmc[{chs}]/O'.format(chs=self.chs))
        # self.outTree.Branch('diaChSeedNoCmc', self.dut_ADC_is_seed, 'diaChSeedNoCmc[{chs}]/O'.format(chs=self.chs))
        # self.outTree.Branch('diaChSignalNoCmc', self.dut_ADC_signal, 'diaChSignalNoCmc[{chs}]/F'.format(chs=self.chs))
        # self.outTree.Branch('diaChPedMeanNoCmc', self.dut_ADC_mean, 'diaChPedMeanNoCmc[{f}]/F'.format(f=self.chs))
        # self.outTree.Branch('diaChPedSigmaNoCmc', self.dut_ADC_sigma, 'diaChPedSigmaNoCmc[{f}]/F'.format(f=self.chs))
        # self.outTree.Branch('diaChCm', self.dut_ADC_cm_chs, 'diaChCm[{f}]/O'.format(f=self.chs))
        self.outTree.Branch('cmn', self.dut_ADC_cm, 'cmn/F')
        # self.outTree.Branch('diaChPed', self.dut_ADC_is_ped_cmc, 'diaChPed[{chs}]/O'.format(chs=self.chs))
        # self.outTree.Branch('diaChHits', self.dut_ADC_is_hit_cmc, 'diaChHits[{chs}]/O'.format(chs=self.chs))
        # self.outTree.Branch('diaChSeed', self.dut_ADC_is_seed_cmc, 'diaChSeed[{chs}]/O'.format(chs=self.chs))
        self.outTree.Branch('diaChSignal', self.dut_ADC_signal_cmc, 'diaChSignal[{chs}]/F'.format(chs=self.chs))
        self.outTree.Branch('diaChPedMeanCmc', self.dut_ADC_mean_cmc, 'diaChPedMeanCmc[{f}]/F'.format(f=self.chs))
        self.outTree.Branch('diaChPedSigmaCmc', self.dut_ADC_sigma_cmc, 'diaChPedSigmaCmc[{f}]/F'.format(f=self.chs))
        self.outTree.Branch('slidingLength', self.fill_slide_leng, 'slidingLength/s')
        # self.b6 = self.outTree.Branch('diaPedChs', self.dut_ADC_is_ped, 'diaPedChs[{chs}]/O'.format(chs=self.chs))
        # self.b7 = self.outTree.Branch('diaHitChs', self.dut_ADC_is_hit, 'diaHitChs[{chs}]/O'.format(chs=self.chs))
        # self.b8 = self.outTree.Branch('diaSeedChs', self.dut_ADC_is_seed, 'diaSeedChs[{chs}]/O'.format(chs=self.chs))
        # self.b9 = self.outTree.Branch('diaPedestalMean', self.dut_ADC_mean, 'diaPedestalMean[{f}]/F'.format(f=self.chs))
        # self.b10 = self.outTree.Branch('diaPedestalSigma', self.dut_ADC_sigma, 'diaPedestalSigma[{f}]/F'.format(f=self.chs))
        # self.b11 = self.outTree.Branch('diaCmChs', self.dut_ADC_cm_chs, 'diaCmChs[{f}]/O'.format(f=self.chs))
        # self.b12 = self.outTree.Branch('diaCm', self.dut_ADC_cm, 'diaCm/F')
        # self.b13 = self.outTree.Branch('diaPedChsCmc', self.dut_ADC_is_ped_cmc, 'diaPedChsCmc[{chs}]/O'.format(chs=self.chs))
        # self.b14 = self.outTree.Branch('diaHitChsCmc', self.dut_ADC_is_hit_cmc, 'diaHitChsCmc[{chs}]/O'.format(chs=self.chs))
        # self.b15 = self.outTree.Branch('diaSeedChsCmc', self.dut_ADC_is_seed_cmc, 'diaSeedChsCmc[{chs}]/O'.format(chs=self.chs))
        # self.b16 = self.outTree.Branch('diaPedestalMeanCmc', self.dut_ADC_mean_cmc, 'diaPedestalMeanCmc[{f}]/F'.format(f=self.chs))
        # self.b17 = self.outTree.Branch('diaPedestalSigmaCmc', self.dut_ADC_sigma_cmc, 'diaPedestalSigmaCmc[{f}]/F'.format(f=self.chs))
        if self.showprogress: print 'Done'

    def FillTree(self):
        if self.showprogress:
            print 'Filling tree...'
            nevents = int(self.ev_end - self.ev_ini + 1)
            bar = CreateProgressBarUtils(nevents)
            bar.start()
        for ev in xrange(self.ev_ini, self.ev_end + 1):
            # self.outTree.GetEntry(ev)
            self.LoadArrays(ev)
            self.outTree.Fill()
            if self.showprogress: bar.update(int(ev - self.ev_ini + 1))
        # self.outTree.Write()
        self.outFile.Write()
        self.outFile.Close()
        if self.showprogress: bar.finish()

    def LoadArrays(self, ev):
        np.putmask(self.dut_ADC_all, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_all_mp.get_obj())[:, ev])
        # np.putmask(self.dut_ADC_is_ped, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_ped_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_is_hit, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_hit_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_is_seed, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_seed_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_signal, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_signal_all_mp.get_obj())[:, ev])
        # np.putmask(self.dut_ADC_mean, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_mean_all_mp.get_obj())[:, ev])
        # np.putmask(self.dut_ADC_sigma, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_sigma_all_mp.get_obj())[:, ev])
        # np.putmask(self.dut_ADC_cm_chs, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_chs_cm_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_signal, np.ones(self.chs, '?'), self.dut_ADC_all - self.dut_ADC_mean)
        self.dut_ADC_cm.fill(np.ctypeslib.as_array(self.dut_ADC_cm_all_mp.get_obj())[ev])
        # np.putmask(self.dut_ADC_is_ped_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_ped_cmc_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_is_hit_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_hit_cmc_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_is_seed_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_is_seed_cmc_all_mp.get_obj())[:, ev].astype('?'))
        # np.putmask(self.dut_ADC_signal_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_signal_cmc_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_mean_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_mean_cmc_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_sigma_cmc, np.ones(self.chs, '?'), np.ctypeslib.as_array(self.dut_ADC_sigma_cmc_all_mp.get_obj())[:, ev])
        np.putmask(self.dut_ADC_signal_cmc, np.ones(self.chs, '?'), self.dut_ADC_all - self.dut_ADC_cm - self.dut_ADC_mean_cmc)
        # for branch in [self.b6, self.b7, self.b8, self.b9, self.b10, self.b11, self.b12, self.b13, self.b14, self.b15, self.b16, self.b17]:
        #     branch.Fill()


def main():
    z = PedestalCalculations()


if __name__ == '__main__':
    main()
