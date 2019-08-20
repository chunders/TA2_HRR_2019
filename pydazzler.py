import numpy as np
import re
import pandas as pd
import json
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.integrate as integrate
from scipy.optimize import minimize_scalar
def omega2lam(w,S_w,axis=0,lAxis=None):
    L = 2*np.pi*299.792/w     
    S_wSize = np.shape(S_w)

    #    S_w = np.expand_dims(S_w,axis=1)
    #    S_wSize = np.shape(S_w)
    if lAxis is None:  
        lAxis = np.linspace(L.min(),L.max(),np.size(L))
    
    if np.size(S_wSize)==2:
        if axis ==0:
            w = np.expand_dims(w,axis=1)
            w2 = np.repeat(w**2, S_wSize[1],axis=1)
        else:
            w = np.expand_dims(w,axis=0)
            w2 = np.repeat(w**2, S_wSize[0],axis=0)
    else:
        w2 = w**2   

    #f_l = interpolate.interp1d(L, S_w*w2/(2*np.pi*299.792),axis=axis)
    f_l = interpolate.interp1d(L, S_w,axis=axis)
    S_l = f_l(lAxis)

    return lAxis, S_l

def lam2omega(l,S_l,axis=0,wAxis=None):
    w = 2*np.pi*299.792/l     
    S_lSize = np.shape(S_l)
    
    #    S_l = np.expand_dims(S_l,axis=1)
    #   S_lSize = np.shape(S_l)
    if wAxis is None:
        wAxis = np.linspace(w.min(),w.max(),np.size(w))
    if np.size(S_lSize)==2:
        if axis ==0:
            l = np.expand_dims(l,axis=1)
            l2 = np.repeat(l**2, S_lSize[1],axis=1)
        else:
            l = np.expand_dims(l,axis=0)
            l2 = np.repeat(l**2, S_lSize[0],axis=0)
    else:
        l2 = l**2
        
    #f_w = interpolate.interp1d(w, S_l*l2/(2*np.pi*299.792),axis=axis)
    f_w = interpolate.interp1d(w, S_l,axis=axis)
    S_w = f_w(wAxis)
    return wAxis, S_w

def E_w2E_t(E_w):
    E_t = (np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(E_w))))
    return E_t
def E_t2E_w(E_t):
    E_w = (np.fft.fftshift(np.fft.fft(np.fft.ifftshift(E_t))))
    return E_w

def normalise(I):
    I = I /np.max(np.abs(I))
    return I


class dazzler( object ): # STYLE: (object) hasn't been needed since Python 3.0
    def __init__( self ):
        self.dazSettings = {}
        self.lam  = []
        self.phi_lam = []
        self.N=2**16
        self.xtalWidth = 8500
        self.phaseLoad=0
    

    def loadFile(self,filepath): # shouldn't this be in __init__?
        """Loads settings from a Dazzler wave file.
        
        Arguments:
         filepath: a string giving a path to a Dazzler wave file
        """

        dazSettings = {}
        lam  = []
        phi_lam = []
        fid = open(filepath) # STYLE: use a with or a try finally
        lCont=1 # STYLE: this could be True
        phaseLoad=0 # STYLE: and this could be False

        while lCont==1: # STYLE: this could be "while lCont:" (even if lCont is an int)
            line = next(fid)
            if line.rfind('!The following is ONLY a reminder of your settings when the wave was saved.')>=0:
                # ^STYLE^: try "if '!The ... was saved.' in line:" instead
                lCont=0 # STYLE: you can use False here, or break explicitly
            elif line.rfind('=')>0: # STYLE: why rfind and not find?
            
                command, value = line.strip().split('=', 1)
                dazSettings[command] = np.float64(value.strip())
            elif line.rfind('#phase')>=0:
                
                phaseLoad =1
            elif phaseLoad==1: # STYLE: again, the "==1" is redundant
                if line.rfind('\t')>0:
                    [a,b] =line.strip().split('\t') # STYLE: get rid of the []; I didn't even know Python allowed them there
                    lam.append(np.float64(a))
                    phi_lam.append(np.float64(b))
                else:
                    lCont=0
                
        fid.close()
        if phaseLoad==1:
            lam= np.array(lam)
            phi_lam = np.array(phi_lam)
        else:
            lam = np.linspace(700,900,num=self.N)
            phi_lam = 0 # shouldn't this be np.zeros_like(lam)?
        self.phaseLoad=phaseLoad    
        self.dazSettings = dazSettings
        self.lam  = lam
        self.phi_lam = phi_lam

    def showParameters(self):
        print(json.dumps(self.dazSettings, indent=2, sort_keys=False))

    def calcWave(self):
        dazSettings = self.dazSettings
        lam = self.lam
        phi_lam = self.phi_lam
        N = self.N
        xtalWidth = self.xtalWidth
        c = 299.792

        # central angular frequency
        w_0 = 2*np.pi*c/(dazSettings['centralwl'])

        # Define the complex function of omega which has the amplitude function  of the dazzler and the phase terms
        w,Amp_w = self.spectralAmp()
        phi_w = self.gatherPhi(w,w_0,dazSettings)

        lambdaTemp, S_lambdaTemp = omega2lam(w,Amp_w)
        S_lambda_fun = interpolate.interp1d(lambdaTemp,S_lambdaTemp,kind='linear')
        S_lambda = S_lambda_fun(lam)
        S_lambda_norm = 1./S_lambda.max()
        S_lambda = S_lambda*S_lambda_norm
        # add phase
        E_w = np.sqrt(Amp_w)*np.exp(-1j*phi_w)
        E_t = E_w2E_t(E_w)
        # time axis
        dt = 2*np.pi/(w.max()-w.min())
        t = np.linspace(-N/2*dt,(N/2-1)*dt,num=N)
                    
       
        E_w_filt = self.applyCrystalClip(w,E_w)
        S_w_filt = np.abs(E_w_filt**2)

        lambdaTemp, S_lambdaTemp_filt = omega2lam(w,S_w_filt)
        S_lambda_filt_fun = interpolate.interp1d(lambdaTemp,S_lambdaTemp_filt,kind='linear')
        S_lambda_filt = S_lambda_filt_fun(lam)*S_lambda_norm           
        
        return lam,S_lambda_filt,S_lambda,t,E_t

    def calcTransErr(self,l,S_l):
        (lam,S_lambda_filt,S_lambda,t,E_t) = self.calcWave()
        S_l_fun = interpolate.interp1d(l,S_l,kind='linear',bounds_error=False, fill_value=0)
        S_l_pulse = S_l_fun(lam)
        #S_trans = np.sum(S_l_pulse*S_lambda_filt)/(np.sum(S_l_pulse*S_lambda))
        S_transErr = np.sqrt(np.sum(((S_lambda_filt-S_lambda)**2)*S_l_pulse)/np.sum(S_l_pulse))
        return S_transErr

    def optimiseDelay(self,l,S_l):
        def rejFun(delay):
            self.dazSettings['delay'] = delay
            return self.calcTransErr(l,S_l)
        res = minimize_scalar(rejFun, [0, self.xtalWidth], options={'xtol': 0.02, 'disp': False})
        self.dazSettings['delay'] = res.x
        return res.x, res.fun

    def plotWaves(self,l=None,S_l=None):
        xtalWidth=self.xtalWidth
        # set up figure
        #fh = plt.figure()
        fh, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
        (lam,S_lambda_filt,S_lambda,t,E_t) = self.calcWave()
        ax1.plot(t,normalise(np.real(E_t)))
        ax1.plot((0, 0),(-1,1),'--k')
        ax1.plot(np.array((1,1))*xtalWidth,(-1,1),'--k')
        ax1.set_xlim(xtalWidth*np.array((-0.5,1.5)))        
        ax1.set_xlabel('Time [ fs ]')
        ax1.set_ylabel('Dazzler wave')
                
        ax2.plot(lam,S_lambda_filt,lam,S_lambda)
        if l is not None:
            ax2.plot(l,S_l)
        ax2.set_ylim((0,1.0))
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_xlabel('Wavelength [ nm ]')
        ax2.set_ylabel('Diffracted amplitude')
        ax2.set_xlim((650, 950))
    
    def calcOpticalPulse(self,PulseObj):
        c = 299.792
        N = PulseObj.N
        dazSettings = self.dazSettings
        pulseSettings = PulseObj.dazSettings
        w = PulseObj.w
        daz_w_0 = 2*np.pi*c/(dazSettings['centralwl'])
        pulse_w_0 = 2*np.pi*c/(pulseSettings['centralwl'])

        # time axis
        dt = 2*np.pi/(w.max()-w.min())
        t = np.linspace(-N/2*dt,(N/2-1)*dt,num=N)
        dwdt = np.mean(np.diff(w))/np.mean(np.diff(t))

        # measured pulse field and phase
        pulse_E_w = PulseObj.E_w
        pulse_phi_W = np.unwrap(np.angle(pulse_E_w))

        # current dazzler phase
        phi_w_daz = self.gatherPhi(w,daz_w_0,dazSettings)
        # dazzler phase for measured pulse
        phi_w_pulse = self.gatherPhi(w,pulse_w_0,PulseObj.dazSettings)

        ## Apply Dazzler Clipping
        # convert to time domain
        S_w = np.abs(pulse_E_w)**2
        S_w = S_w/integrate.trapz(S_w,x=w)*PulseObj.energy
        E_w = np.sqrt(S_w)*np.exp(-1j*(phi_w_daz))
  
        E_w_filt = self.applyCrystalClip(w,E_w)

        # add phase
        E_w_final = np.abs(E_w_filt)*np.exp(-1j*(phi_w_pulse-phi_w_daz))
        E_t_final = E_w2E_t(E_w_final)*np.sqrt(N*dwdt)
        return t, E_t_final

    
    def gatherPhi(self,w,w_0,dazSettings):
        c = 299.792
        phaseLoad = self.phaseLoad
        lam = self.lam
        phi_lam = self.phi_lam
        xtalWidth = self.xtalWidth
        # convert phase to frequency domain
        if phaseLoad==1:
            phi_w_fun = interpolate.interp1d(2*np.pi*c/lam,phi_lam,kind='linear',bounds_error=False, fill_value=0)
            phi_w = phi_w_fun(w)
        else:
            phi_w=0

        # add on 4 orders from dazzler settings
        phi_w = phi_w + dazSettings['delay']*(w-w_0)
        phi_w = phi_w +(dazSettings['order2']/2*(w-w_0)**2)
        phi_w = phi_w +(dazSettings['order3']/6*(w-w_0)**3)
        phi_w = phi_w +(dazSettings['order4']/24*(w-w_0)**4)
        return phi_w

    def applyCrystalClip(self,w,E_w):
        N = np.size(w)
        # convert to time domain
        E_t = E_w2E_t(E_w)
        # time axis
        dt = 2*np.pi/(w.max()-w.min())
        t = np.linspace(-N/2*dt,(N/2-1)*dt,num=N)
                    
        # clip with the temporal window of the dazzler
        E_t_filt = E_t*(t>=0)*(t<=self.xtalWidth)
        # transform back and take the amplitude
        E_w_filt = E_t2E_w(E_t_filt)
        return E_w_filt

    def spectralAmp(self):
        dazSettings = self.dazSettings
        lam = self.lam
        N = self.N
        c = 299.792
        # central angular frequency
        w_0 = 2*np.pi*c/(dazSettings['centralwl'])
        # angular frequency axis
        w = np.linspace(2*np.pi*c/lam.max(),2*np.pi*c/lam.min(),num=N)
        w = np.linspace(0.5,5,num=N)
        phi_w = self.gatherPhi(w,w_0,dazSettings)

        # Define the complex function of omega which has the amplitude function  of the dazzler and the phase terms
        f_w0 = 2*np.pi*c/dazSettings['position']
        xi_0 = dazSettings['width']/(2*dazSettings['position'])
        g_w0 = 2*np.pi*c/dazSettings['hposition']
        xi_0 = dazSettings['width']/(2*dazSettings['position'])
        xi_1 = dazSettings['hwidth']/(2*dazSettings['hposition'])
        f_dw = f_w0*(xi_0-xi_0**3)
        g_dw = g_w0*(xi_1-xi_1**3)/2
        f_w = np.exp(-((w-f_w0)/f_dw)**6)
        with np.errstate(divide='ignore',invalid='ignore'):
            g_w = 1-dazSettings['hdepth']*np.exp(-((w-g_w0)/g_dw)**2)
        Amp_w = f_w*g_w

        return w,Amp_w
    

class laserPulse( object ):
    def __init__( self ):
        self.dazSettings = {}
        self.lam  = []
        self.phi_lam = []
        self.N=2**16
        self.xtalWidth = 8500
        self.phaseLoad=0
        self.energy = 1
        self.w = []
        self.E_w = []
    

    def loadFile(self,filepath):
        dazSettings = {}
        lam  = []
        phi_lam = []
        fid = open(filepath)
        lCont=1
        phaseLoad=0

        while lCont==1:
            line = next(fid)
            if line.rfind('!The following is ONLY a reminder of your settings when the wave was saved.')>=0:
                lCont=0
            elif line.rfind('=')>0:
            
                command, value = line.strip().split('=', 1)
                dazSettings[command] = np.float64(value.strip())
            elif line.rfind('#phase')>=0:
                
                phaseLoad =1
            elif phaseLoad==1:
                if line.rfind('\t')>0:
                    [a,b] =line.strip().split('\t')
                    lam.append(np.float64(a))
                    phi_lam.append(np.float64(b))
                else:
                    lCont=0
                
        fid.close()
        if phaseLoad==1:
            lam= np.array(lam)
            phi_lam = np.array(phi_lam)
        else:
            lam = np.linspace(700,900,num=self.N)
            phi_lam = 0
        self.phaseLoad=phaseLoad    
        self.dazSettings = dazSettings
        self.lam  = lam
        self.phi_lam = phi_lam


    def showParameters(self):
        print(json.dumps(self.dazSettings, indent=2, sort_keys=False))

    