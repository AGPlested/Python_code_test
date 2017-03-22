
#import pyhht as HHT

import matplotlib as MP
import matplotlib.pyplot as MPP
import pylab as PL

import re as RE
import pandas as PD
import numpy as NP
import scipy.signal as SIG
#import scipy.ndimage as NDI

from lmfit import minimize as FIT, Parameters as FIT_Params
# import openpyxl as XL

import axographio as AG
import os as OS


### Coverage:
"""
    - dictionaries
    - lambda
    - Objects
    - axograph io
    - data frames
    - plotting & LaTeX
    - fiting

"""
### possible ToDos:
"""
    - fit double gaussian
    - empirical mode decomposition
    - plot smoothed curve (NDI.gaussian_filter1d)

"""


###########################
###   Global functions  ###
###    and definitions  ###
###########################

def ParamsToDict(params):
    params_dict = {}
    for key in params.keys():
        params_dict[key] = params[key].value
    return params_dict



rmse = lambda params, fcn, x, data:  \
                    NP.abs(data - fcn(ParamsToDict(params), x))


units = {   '1': 1 \
            , 'ms': 1e3 \
            , 'pA': 1e12 \
            }


the_font = {'family':'sans-serif'
            ,'sans-serif':['Liberation Sans'] # Helvetica # Liberation Sans
            ,'weight' : 'normal' \
            ,'size'   : 12 \
            }



###########################
###   Trace Object      ###
###########################

class AxoTrace(object):
    # CONSTRUCTOR:
    def __init__(self, file_path, trace_nr = 0):
    # load trace from axograph file
        self.file = file_path

        self.date = RE.search(r'([0-9]{6})', self.file ).group(0)
        self.filenum = int(RE.search(r' ([0-9]{3}) ', self.file ).group(0))
        self.episode = RE.sub('ep', 'episode ', RE.search(r'ep([0-9]{2})', self.file ).group(0) )


        axodata = AG.read(self.file)
    # extract time
        times = [not(elmt is None) for elmt in [RE.search(r'(?i)(Time)', nstr) for nstr in axodata.names]]
        trace_t = (NP.array([axodata.data[idx2][:] for idx2 in [idx for idx in range(len(times)) if times[idx]] ] , dtype = float)[0]) * units['ms']

    # extract current
        restr = r'(?i)(Ipatch)' # current
        da = [not(elmt is None) for elmt in [RE.search(restr, nstr) \
                            for nstr in axodata.names]]
        trace_i = [NP.array(axodata.data[idx2]) \
                            for idx2 in [idx for idx in range(len(da)) if da[idx]] ]

        pd_dict = { \
                      "time": trace_t \
                    , "current": NP.array(trace_i[trace_nr]) * units['pA'] \
                    }

        self.raw = PD.DataFrame.from_dict(pd_dict)
        # print self.raw[self.raw.time <= 1].current

        self.sampling_rate = (len(self.raw.current)-1) / ((NP.max(self.raw.time)-NP.min(self.raw.time))/units['ms'])
        print self.sampling_rate 
            
            #+trace_t[:-1]/2
        pd_dict = { \
                      "time": trace_t[1:] \
                    , "current": NP.diff(self.raw.current) \
                    }
        self.diff = PD.DataFrame.from_dict(pd_dict)


    def PlotTrace(self, ax):
        ax.plot(self.raw.time, self.raw.current, 'k-')

        ax.set_xlabel(r'time $\left( ms\right)$')
        ax.set_ylabel(r'current $\left( pA\right)$')


    def PlotDiffTrace(self, ax):

        ax.plot(self.diff.time, self.diff.current, 'k-')

        ax.set_title(r'gradient $\nabla \mathrm{I}$')
        ax.set_xlabel(r'time $\left( ms\right)$')
        ax.set_ylabel(r'$\Delta$I $\left( pA\right)$')


    def HistoDiff(self, ax):
        
        data = NP.random.normal(size=10000)
        PL.hist(data, bins=NP.logspace(0.1, 1.0, 50))
        PL.gca().set_xscale("log")
        PL.show()

        #hist, edges = NP.histogram(a = (self.raw.current), bins = 50, density = True)
        #width = 0.9 * (edges[1] - edges[0])
        #bin_centers = (edges[1:]+edges[:-1])/2
        #ax.bar(bin_centers, hist, align='center'\
         #               , width=width, edgecolor = (0,0,0), facecolor = (1,1,1)\
          #              , label = r'binned data'\
           #             )


        # fit with single gauss
        params = FIT_Params()
        params.add('a1', value=0.8, min = 0, max = 1) # , min = 0
        params.add('mu1', value=0, min = -5000, max = 5000) # , min = 0
        params.add('sigma1', value=1, min = 1e-6, max = 5000) # , min = 0
        params.add('a2', value=0.2, min = 0, max = 1) # , min = 0
        params.add('mu2', value=-5, min = -5000, max = 5000) # , min = 0
        params.add('sigma2', value=1, min = 1e-6, max = 5000) # , min = 0

        # print ParamsToDict(params)
        gauss = lambda p, x: \
                        p['a1']*NP.exp(-((x-p['mu1'])**2)/(2*(p['sigma1']**2)))/(p['sigma1']*NP.sqrt(2*NP.pi)) \
                        + p['a2']*NP.exp(-((x-p['mu2'])**2)/(2*(p['sigma2']**2)))/(p['sigma2']*NP.sqrt(2*NP.pi))

        # TODO: fit with double gaussian!


        fit_fcn = gauss

        out = FIT(rmse, params, args=( fit_fcn \
                                    , bin_centers \
                                    , hist \
                                    ) \
                                , method = 'leastsq' \
                                #, maxfev= 10000    # optional for "leastsq"
                                #, tol = 1e-12      # optional for "nelder"
                                ) 

        print "fit results: %s" % out.values
        ax.plot(bin_centers, gauss(out.values, bin_centers), 'r:', label = r'mono-gaussian fit')

        # polish plot
        ax.set_title(r'histogram of $\nabla \mathrm{I}$')
        ax.set_xlabel(r'current step')
        ax.set_ylabel(r'probability')

        ax.legend(loc=1)



    def EMD(self, cycles = 3):
        pass
        # TODO: implement this



###########################
###   Main process      ###
###########################
def main():
### Load data ###
    trace = AxoTrace( \
                        "140905 001 Average-ep61-2kHz-filt.axgx" \
                        # "/Users/Jelena/files/140905 001 Average-ep61-not-filt.axgx" \
                    ) 
    # print trace.sampling_rate

    

### HHT ###
    # TODO: decompose to IMF




### figure ###

    figwidth = 16/2.54*1.25
    figheight = 16/2.54*1.25

# define figure
    fig = MPP.figure(figsize=(figwidth, figheight), facecolor = None, dpi = 150)
    fig.hold(True) #(if multiple figures)
    fig.subplots_adjust(top = 0.95, right=0.95, bottom = 0.1, left = 0.06, wspace=0.2, hspace=0.2)

    print ( trace.filenum )
    fig.suptitle( "%s, file %04.0f, %s" % ( trace.date, trace.filenum, trace.episode ) )

# define subplots
    gs = MP.gridspec.GridSpec(2, 2) #  height_ratios=[2, 3]
    axarr = []
    axarr.append(fig.add_subplot(gs[0, :]))
    axarr.append(fig.add_subplot(gs[1, 0], sharex=axarr[0]))
    axarr.append(fig.add_subplot(gs[1, 1]))

    for ax in axarr:
            ax.get_xaxis().set_tick_params(which='both', direction='out')
            ax.get_yaxis().set_tick_params(which='both', direction='out')
            ax.tick_params(top="off")
            ax.tick_params(right="off")
            
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(False)

            # # log scaled axes:
            #ax.set_xscale('log')
            #ax.xaxis.set_major_formatter(MP.ticker.ScalarFormatter())


# plot trace etc.
    trace.PlotTrace(axarr[0])
    trace.PlotDiffTrace(axarr[1])
    trace.HistoDiff(axarr[2])


#LaTeX in python plots:
    MP.rc( 'font', **the_font)

    MPP.rcParams.update(**{'text.usetex': True})

    MP.rcParams['text.latex.preamble'] = [\
            r'\usepackage{upgreek}'
            , r'\usepackage{cmbright}'
            , r'\usepackage{sansmath}' 
            , r'\sansmath' 
            ]
    
# save as svg (or pdf)
    fig.savefig( RE.sub('.axgx', '.svg',  "%s/%s" % ( OS.getcwd(), trace.file  ) )\
                            , transparent=False, dpi=150, format = 'svg')


# don't forget to show!
    MPP.show()



if __name__ == "__main__":
    # this is here because all the function definitions should happen before processing the main algorithm.
    main()


