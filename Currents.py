#!/usr/bin/env python
# --------------------------------------------------------
#       Script to show the recorded current for the 3D beam test runs
# created on September 4th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from glob import glob
import time
import ipdb
import cPickle as pickle

import json as js
from collections import OrderedDict
import os, sys, shutil

from ConfigParser import ConfigParser
from optparse import OptionParser

import numpy as np
import ROOT as ro
from Utils import *

# ====================================
# CONSTANTS
# ====================================
evts_per_file = 10000
tree_name = 'hv_currents'
# ====================================
# CLASS FOR THE DATA
# ====================================
class Currents:
    """reads in information from the keithley log file"""

    def __init__(self, settings_f, diam='', pola='both'):

        # settings
        self.settings_file = Correct_Path(settings_f)
        self.diamond = diam
        self.polarity = pola
        self.time_offset = 0.0

        # config
        self.currents_logs_dir = ''
        self.testbeam_log_file = ''
        self.testbeam_log_file_name = ''
        self.out_dir = ''
        self.subdir = 'Currents'
        self.raw_files_dir = ''
        self.ReadSettingsFile()

        self.runs_evets_times = {}
        self.runs_info = None
        self.runs = []
        self.LoadRunsInfo()

        self.hv_channels = []
        self.hv_channels_root = ro.TFile()
        self.hv_channels_tree = ro.TTree()
        self.CheckHVChannelsROOTFiles()

        self.Histos = {}
        self.Stuff = []


    def ReadSettingsFile(self):
        if os.path.isfile(self.settings_file):
            pars = ConfigParser()
            pars.read(self.settings_file)
            print 'Loading settings file:', self.settings_file

            if pars.has_section('SETTINGS'):
                if pars.has_option('SETTINGS', 'currents_logs_dir'):
                    self.currents_logs_dir = pars.get('SETTINGS', 'currents_logs_dir')
                else:
                    print 'The option "currents_logs_dir" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'testbeam_log_file'):
                    self.testbeam_log_file = pars.get('SETTINGS', 'testbeam_log_file')
                    self.testbeam_log_file_name = self.testbeam_log_file.split('/')[-1]
                else:
                    print 'The option "testbeam_log_file" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'raw_files_dir'):
                    self.raw_files_dir = pars.get('SETTINGS', 'raw_files_dir')
                else:
                    print 'The option "raw_files_dir" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'outputdir'):
                    self.out_dir = pars.get('SETTINGS', 'outputdir')
                else:
                    print 'The option "outputdir" under the section "SETTINGS" is required'
                if pars.has_option('SETTINGS', 'subdir'):
                    self.subdir = pars.get('SETTINGS', 'subdir')
            else:
                print 'The section "SETTINGS" is needed in the settings file'

    # ==========================================================================
    # region INIT

    def LoadRunsInfo(self, force_recreate=False):
        self.runs_info = self.GetRunInfoJson()
        self.runs_evets_times = self.GetEventsTimesPickle()
        if (not self.runs_info) or (force_recreate):
            if self.testbeam_log_file != '':
                temp_runs_info = OrderedDict()
                with open(self.testbeam_log_file, 'r') as ftb:
                    lines = ftb.readlines()
                    for line in lines:
                        if len(line.split()) > 1 and (IsInt(line.split()[0]) or ('/' in line.split()[0] and IsInt(line.split()[1]))):
                            run_i = str(line.split()[0]) if IsInt(line.split()[0]) else str(line.split()[1])
                            diam_i = str(line.split()[1]) if IsInt(line.split()[0]) else str(line.split()[2])
                            volt_i = line.split()[2].strip('V') if IsInt(line.split()[0]) else line.split()[3].strip('V')
                            if '-' in volt_i:
                                volt_i = volt_i.strip('-')
                                volt_i = -1000 * int(volt_i.strip('K').strip('k')) if 'K' in volt_i or 'k' in volt_i else -1 * int(volt_i)
                            else:
                                volt_i = volt_i.strip('+')
                                volt_i = 1000 * int(volt_i.strip('K').strip('k')) if 'K' in volt_i or 'k' in volt_i else int(volt_i)
                            numev_i = line.split()[8] if IsInt(line.split()[0]) else line.split()[9]
                            numev_i = 1000 * int(numev_i.strip('K').strip('k')) if 'K' in numev_i or 'k' in numev_i else int(numev_i)
                            hv_folder_i = Get_From_User_Line('HV folder as it is in the currents logs directory for run ' + run_i + ' with voltage ' + str(volt_i) + ' and DUT ' + diam_i, update=False)
                            temp_runs_info[run_i] = {'dia': diam_i, 'bias': volt_i, 'num_ev': numev_i, 't_ini': -1.0, 't_end': -1.0, 'hv_ch': hv_folder_i}
                            print run_i, diam_i, volt_i, numev_i, hv_folder_i
                self.runs_info = temp_runs_info
                self.GetDictOfRawFiles()
                self.DumpRunsInfoInJson(force_recreate)
                # with open('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name), 'w') as fj:
                #     js.dump(temp_runs_info, fj, indent=2, sort_keys=True)
                #     fj.truncate()
                # print '{c} runs info and saved it in json: {d}/{f}.json'.format(c='Created' if not self.runs_info else 'Recreated', d=self.currents_logs_dir, f=self.testbeam_log_file_name)
        elif self.runs_info:
            print 'Loaded runs info from json: {d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name)

    def GetRunInfoJson(self):
        if os.path.isfile('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name)):
            with open('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name), 'r') as fj:
                temp_runs_info = js.load(fj)
            return OrderedDict(sorted(temp_runs_info.iteritems()))

    def GetEventsTimesPickle(self):
        temp = {}
        if os.path.isfile('{d}/{f}.pkl'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name)):
            with open('{d}/{f}.pkl'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name), 'rb') as pkl:
                temp = pickle.load(pkl)
                print 'Loaded pickle file with the event times for all the runs in the testbeam'
        return temp


    def RecreateRunsInfo(self):
        self.LoadRunsInfo(True)

    def DumpRunsInfoInJson(self, isrecreated=False):
        with open('{d}/{f}.json'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name), 'w') as fj:
            js.dump(self.runs_info, fj, indent=2, sort_keys=True)
            fj.truncate()
        print '{c} runs info and saved it in json: {d}/{f}.json'.format(c='Created' if not isrecreated else 'Recreated', d=self.currents_logs_dir, f=self.testbeam_log_file_name)

    def GetDictOfRawFiles(self):
        if os.path.isdir(self.raw_files_dir):
            runs_paths = glob(self.raw_files_dir + '/*')
            self.runs = sorted([run_i.split('/')[-1] for run_i in runs_paths if IsInt(run_i.split('/')[-1])])
            for run_i in self.runs:
                if run_i in self.runs_info.keys():
                    raws_run_i = [raw_i for raw_i in glob('{d}/{r}/RUN*'.format(d=self.raw_files_dir, r=run_i)) if IsInt(raw_i.split('.')[-2].split('_')[-1])]
                    raws_run_i.sort(key=lambda x: os.path.getmtime(x))
                    raws_run_i_tf = [os.path.getmtime(raw_i) + self.time_offset for raw_i in raws_run_i]
                    raws_run_i_ti = [raws_run_i_tf[i - 1] + self.time_offset for i in xrange(1, len(raws_run_i))]
                    try:
                        raws_run_i_ti.insert(0, float(RoundInt((len(raws_run_i) * raws_run_i_tf[0] - raws_run_i_tf[-1]) / float(len(raws_run_i) - 1))))
                    except Exception:
                        print run_i
                        ipdb.set_trace()
                    self.runs_evets_times[run_i] = OrderedDict()
                    # for pos_i, [ti_i, tf_i] in enumerate(zip(raws_run_i_ti, raws_run_i_tf)):
                    #     if (pos_i + 1) * evts_per_file < self.runs_info[run_i]['num_ev']:
                    #         for ev in xrange(pos_i * evts_per_file, (pos_i + 1) * evts_per_file):
                    #             evts_times[ev] = (evts_per_file * ti_i + (ev % evts_per_file) * (tf_i - ti_i)) / float(evts_per_file)
                    #     else:
                    #         delta_t_last = (tf_i - ti_i) / float((self.runs_info[run_i]['num_ev'] - 1) - pos_i * evts_per_file)
                    #         for ev in xrange(pos_i * evts_per_file, self.runs_info[run_i]['num_ev']):
                    #             evts_times[ev] = ti_i + (ev % evts_per_file) * delta_t_last
                    self.runs_info[run_i]['t_ini'] = raws_run_i_ti[0]
                    self.runs_info[run_i]['t_end'] = raws_run_i_tf[-1]
                    for pos_i, [ti_i, tf_i] in enumerate(zip(raws_run_i_ti, raws_run_i_tf)):
                        delta_t = (tf_i - ti_i) / float(evts_per_file) if (pos_i + 1) * evts_per_file < self.runs_info[run_i]['num_ev'] else (tf_i - ti_i) / float((self.runs_info[run_i]['num_ev'] - 1) - pos_i * evts_per_file)
                        lim_ev_sup = (pos_i + 1) * evts_per_file if (pos_i + 1) * evts_per_file < self.runs_info[run_i]['num_ev'] else self.runs_info[run_i]['num_ev']
                        for ev in xrange(pos_i * evts_per_file, lim_ev_sup):
                            self.runs_evets_times[run_i][ev] = ti_i + (ev % evts_per_file) * delta_t

            self.DumpPickleRunsEventTimes()

    def DumpPickleRunsEventTimes(self):
        if len(self.runs_evets_times) > 0:
            picklepath = '{d}/{f}.pkl'.format(d=self.currents_logs_dir, f=self.testbeam_log_file_name)
            pickle.dump(self.runs_evets_times, open(picklepath, 'wb'))
            print 'Created pickle file with the event times for all the runs in the testbeam'


    def CheckHVChannelsROOTFiles(self, load_option='READ'):
        self.hv_channels = sorted([dir_i for dir_i in os.listdir(self.currents_logs_dir) if dir_i in list(set([elem_i['hv_ch'] for elem_i in self.runs_info.values()]))])
        for hv_ch in self.hv_channels:
            if load_option.lower() in ['recreate', 'new']:
                self.hv_channels_tree = None
                self.hv_channels_root = ro.TFile(self.currents_logs_dir + '/' + hv_ch + '.root', load_option)
                self.hv_channels_tree = ro.TTree(tree_name, tree_name)
                self.FillCurrentsTreeOfHVChannelFolder(hv_ch)

    def FillCurrentsTreeOfHVChannelFolder(self, hvfolder):
        print 'Filling {t} tree for HV folder {f}'.format(t=self.hv_channels_tree, f=hvfolder)
        if os.path.isdir(self.currents_logs_dir + '/' + hvfolder):
            files = glob(self.currents_logs_dir + '/' + hvfolder + '/*.log')
            files.sort(key=lambda x: os.path.getmtime(x))
            # create branches
            hv_ch = ro.TString(hvfolder)
            # hv_ch_time = ro.TDatime()
            hv_ch_time = ro.TTimeStamp()
            hv_ch_volt = np.zeros(1, 'f4')
            hv_ch_curr = np.zeros(1, 'f4')
            hv_ch_run = np.zeros(1, 'uint16')
            # hv_ch_ev_i = np.zeros(1, 'int32')
            # hv_ch_ev_f = np.zeros(1, 'int32')
            hv_ch_dia = ro.TString()
            # assign branches to tree
            self.hv_channels_tree.Branch('channelHV', hv_ch)
            self.hv_channels_tree.Branch('voltageHV', hv_ch_volt, 'voltageHV/F')
            self.hv_channels_tree.Branch('currentHV', hv_ch_curr, 'currentHV/F')
            self.hv_channels_tree.Branch('timeHV', hv_ch_time)
            self.hv_channels_tree.Branch('runHV', hv_ch_run, 'runHV/s')
            # self.hv_channels_tree.Branch('eventIni', hv_ch_ev_i, 'eventIni/I')
            # self.hv_channels_tree.Branch('eventFin', hv_ch_ev_f, 'eventFin/I')
            self.hv_channels_tree.Branch('diaHV', hv_ch_dia)

            print 'Reading all files inside', hvfolder
            lines_f = []
            bar = CreateProgressBarUtils(len(files))
            bar.start()
            for it, file_i in enumerate(files):
                tempf = file_i.strip('.log').split('_')
                year_f, month_f, day_f = tempf[-6], tempf[-5], tempf[-4]
                with open(file_i, 'r') as logf:
                    for line in logf:
                        lines_f.append('{y} {m} {d} {l}'.format(y=year_f, m=month_f, d=day_f, l=line))
                    bar.update(it + 1)
            bar.finish()
            print 'Removing lines not used...', ; sys.stdout.flush()
            list_f = [[int(line.split()[0]), int(line.split()[1]), int(line.split()[2]), line.split()[3], float(line.split()[4]), float(line.split()[5])] for line in lines_f if len(line.strip('\n').split()) > 5 and IsFloat(line.split()[4]) and IsFloat(line.split()[5])]
            print 'Done'
            print 'Filling tree:'
            bar = CreateProgressBarUtils(len(list_f))
            bar.start()
            for it, line in enumerate(list_f):
                hv_ch_dia.Clear()
                year_f, month_f, day_f = int(line[0]), int(line[1]), int(line[2])
                HMS = line[3].split(':')
                hour_f, min_f, sec_f = int(HMS[0]), int(HMS[1]), int(HMS[2])
                volt_f, curr_f = line[4], line[5]
                # hv_ch_time.Set(year_f, month_f, day_f, hour_f, min_f, sec_f)
                hv_ch_time.Set(year_f, month_f, day_f, hour_f, min_f, sec_f, 0, False, 0)
                run_f = self.GuessRunFromRunInfo(hvfolder, hv_ch_time.AsDouble())
                dia_f = self.runs_info[run_f]['dia'] if run_f in self.runs_info.keys() else 'None'
                # evi_f, evf_f = self.FindEventIniAndFinFromTimeStamp(hv_ch_time.Convert(), run_f)
                # evi_f, evf_f = -1, -1
                # evi_f, evf_f = self.FindEventIniAndFinFromTimeStamp(hv_ch_time.AsDouble(), run_f)
                hv_ch_volt.fill(volt_f)
                hv_ch_curr.fill(curr_f)
                hv_ch_run.fill(int(run_f))
                # hv_ch_ev_i.fill(evi_f)
                # hv_ch_ev_f.fill(evf_f)
                hv_ch_dia.Append(dia_f)
                self.hv_channels_tree.Fill()
                bar.update(it + 1)

            #
            # for it, file_i in enumerate(files):
            #     with open(file_i, 'r') as logf:
            #         lines_f = logf.readlines()
            #         # remove lines not used
            #         list_f = [[line.strip('\n').split()[0], float(line.strip('\n').split()[1]), float(line.strip('\n').split()[2])] for line in lines_f if len(line.strip('\n').split()) > 2 and IsFloat(line.strip('\n').split()[1]) and IsFloat(line.strip('\n').split()[2])]
            #         hv_ch_dia.Clear()
            #         for line in list_f:
            #             HMS = line[0].split(':')
            #             hour_f, min_f, sec_f = int(HMS[0]), int(HMS[1]), int(HMS[2])
            #             volt_f, curr_f = line[1], line[2]
            #             hv_ch_time.Set(year_f, month_f, day_f, hour_f, min_f, sec_f, 0, False, 0)
            #             run_f = self.GuessRunFromRunInfo(hvfolder, hv_ch_time)
            #             dia_f = self.runs_info[str(run_f)]['dia'] if str(run_f) in self.runs_info.keys() else 'None'
            #             evi_f, evf_f = self.FindEventIniAndFinFromTimeStamp(hv_ch_time, run_f)
            #             hv_ch_volt.fill(volt_f)
            #             hv_ch_curr.fill(curr_f)
            #             hv_ch_run.fill(run_f)
            #             hv_ch_ev_i.fill(evi_f)
            #             hv_ch_ev_f.fill(evf_f)
            #             hv_ch_dia.Append(dia_f)
            #             self.hv_channels_tree.Fill()
            #             bar.update(it + 1)
            self.hv_channels_root.Write()
            self.hv_channels_root.Close()
            bar.finish()

    def GuessRunFromRunInfo(self, hvfolder, timestamp):
        for run_i, items in self.runs_info.iteritems():
            if hvfolder in items['hv_ch']:
                if items['t_ini'] <= timestamp <= items['t_end']:
                    return run_i
        return '0'

    def FindEventIniAndFinFromTimeStamp(self, timestamp, run_f='0'):
        if run_f in self.runs_info.keys():
            cond = np.less_equal(np.abs(np.subtract(timestamp, self.runs_evets_times[run_f].values())), 0.5)
            return cond.argmax(), self.runs_info[run_f]['num_ev'] - 1 - cond[::-1].argmax()
        return -1, -1

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings_f', type='string', help='Settings file containing the information of the currents and hv runlogs')
    parser.add_option('-d', '--diamond', dest='diamond', type='string', default='', help='diamond to be analysed. If not given, the default is empty string')
    parser.add_option('-p', '--polarity', dest='polarity', type='string', default='both', help='selects the polarity to analyse. The valid options are: negative, positive, both. The default is "both"')
    parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true', help='enables automatic analysis with a given diamond')

    (options, args) = parser.parse_args()
    settings_f = str(options.settings_f)
    diam = str(options.diamond)
    pola = str(options.polarity)
    autom = bool(options.automatic)

    cu = Currents(settings_f, diam, pola)

    if autom:
        pass
