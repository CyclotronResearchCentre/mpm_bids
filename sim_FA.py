import numpy as np
import configparser
import re
import nibabel as nib
import matplotlib.pyplot as plt

def Rz(phi):
    return np.matrix([[np.cos(phi), -np.sin(phi), 0 ],[np.sin(phi), np.cos(phi),0],[0,0,1]])

def Rp(Omega,omega1,t):
    omegae = np.sqrt(Omega*Omega+omega1*omega1)
    cw = np.cos(omegae*t); sw = np.sin(omegae*t)
    cp = omega1/omegae; sp = Omega/omegae
    return np.matrix([[cp*cp+cw*sp*sp,   sp*sw,  cp*sp*(1-cw)],
                      [-sw*sp,           cw,     sw*cp],
                      [cp*sp*(1-cw),     -sw*cp, sp*sp+cw*cp*cp]])
    
def R(phi,Omega,omega1,t):
    if omega1==0:
        omega1=1E-12
    return Rz(phi)@Rp(Omega,omega1,t)@Rz(-phi)

class Pulse_from_ini:

    def __init__(self,iniFile,B1map):
        self.iniFile = iniFile
        self.B1map   = B1map

        self.read_B1map()
        self.read_ini()
        self.calc_Pulse()


    def read_B1map(self):
        # strip ".nii"  from filename
        f = self.B1map[:-4]
        
        # load data
        mag = nib.load(f+".nii")
        mag_data = mag.get_fdata().astype(np.float32)

        pha = nib.load(f+"_ph.nii")
        pha_data = pha.get_fdata().astype(np.float32)

        self.B1 = mag_data*np.exp(1j*(pha_data-2048)/1800*np.pi)

        self.affine = mag.affine
        self.header = mag.header

        print(self.B1.shape, self.B1.dtype)

    def read_ini(self):
        config = configparser.ConfigParser()
        config.read(self.iniFile)

        self.deltaT = int(re.findall(r'\b\d+\b', config["Gradient"]["GradRasterTime"])[0])
        self.points = int(re.findall(r'\b\d+\b', config["Gradient"]["GradientSamples"])[0])
        self.nomFA  = int(re.findall(r'\b\d+\b', config["pTXPulse"]["NominalFlipAngle"])[0])

        print(self.deltaT)

        self.Volt_mag = np.zeros((self.points,8))
        self.Volt_pha = np.zeros((self.points,8))
        self.Grad     = np.zeros((self.points,3))

        for c in range(8):
            for i in range(self.points):
                self.Volt_mag[i,c], self.Volt_pha[i,c] = config["pTXPulse_ch%i"%c]["rf[%i]"%i].split()
                #self.Volt_pha[i,c] = np.pi/4*c

        for i in range(self.points):
            self.Grad[i,1],self.Grad[i,0],self.Grad[i,2] = config["Gradient"]["G[%i]"%i].split()

    def calc_Pulse(self):
        self.Pulse = self.B1@(self.Volt_mag*np.exp(1j* self.Volt_pha)).transpose()

        self.Pulse *= 4.5366E-5*267E-6
        print(abs(np.max(self.Pulse)))
        ##TODO: rescale

    def calc_omega(self,x,y,z):
        Pos = (self.affine@np.matrix([x,y,z,1]).transpose())[:3]
        Pos[0] *= -1
        O = -np.squeeze(np.array(self.Grad@Pos*267.515)/1E6)
        return O

    def calc_FA(self,x,y,z,offcenter):
        Mz = np.matrix([[0],[0],[1]])

        offset = self.calc_omega(x,y,z)
        
        pulse_local = self.Pulse
        for i in range(self.points):
            j = self.points-i-1
            Rstep = R(np.angle(pulse_local[x,y,z,j]),offset[j]+offcenter,abs(pulse_local[x,y,z,j]),self.deltaT)
            Mz = Rstep@Mz
        return 180*np.arccos(Mz[2])/np.pi, 180*np.arctan2(Mz[0],Mz[1])/np.pi

