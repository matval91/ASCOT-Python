import sys
# GUI 
import Tkinter as tk
import tksty
import tkMessageBox

#retrieve data from tree
import geo,read_sig

#nubeam-related
import nb_wrnml, nbrun

#other stuff
import os,fnmatch,shutil
#import matplotlib.pylab as plt


sty=tksty.TKSTY()

class mkrun:

  def __init__(self, master):
    #if __name__ == '__main__':
    #  self.mainframe = tk.TK()
    #else:
    #  self.mainframe = tk.Toplevel()
    
    self.indict = {}
    #======================
    #FRAME FOR SHOT SETTINGS
    shotframe = sty.mySubframe(master)
    sty.myLabel(shotframe,'shot').grid(row=1, column=1)
    sty.myLabel(shotframe,'dt')  .grid(row=1, column=3)
    sty.myLabel(shotframe,'tbeg').grid(row=2, column=1)
    sty.myLabel(shotframe,'tend').grid(row=2, column=3)
    sty.myLabel(shotframe,'nrho').grid(row=3, column=1)
    sty.myLabel(shotframe,'nthe').grid(row=3, column=3)
    self.eshot = sty.myEntry(shotframe,'53778')
    self.edt   = sty.myEntry(shotframe,'0.02' )
    self.etbeg = sty.myEntry(shotframe,'1.4'   )
    self.etend = sty.myEntry(shotframe,'1.4'  )
    self.enrho = sty.myEntry(shotframe,'41'   )
    self.enthe = sty.myEntry(shotframe,'101'  )
    
    self.eshot.grid(row = 1, column = 2)
    self.edt  .grid(row = 1, column = 4)
    self.etbeg.grid(row = 2, column = 2)
    self.etend.grid(row = 2, column = 4)
    self.enrho.grid(row = 3, column = 2)
    self.enthe.grid(row = 3, column = 4)
    dictbutt = sty.myButton(shotframe,'Dict', self.init_dict, nrow=2, ncol=5)

    #FRAME FOR FUNCTIONS (RETRIEVE DATA)
    commframe = sty.mySubframe(master)
    comm={ 'Namelist':self.write_nml,\
           '2D':self.dw2read, 'GEO':self.geofread, \
          'Run NUBEAM':self.nbrun}
    n_row = 0
    nml_but =sty.myButton(commframe,'Namelist',comm['Namelist'])
    nml_but.grid(row=1, column=1)
    
    dw1 = sty.myButton(commframe,'1D',self.dw1read)
    dw1.grid(row=2, column=1)
    dw2 = sty.myButton(commframe,'2D',self.dw2read)
    dw2.grid(row=2, column=2)
    geo=sty.myButton(commframe,'GEO',comm['GEO'])
    geo.grid(row=2, column=3)

    # FRAME FOR QUIT BUTTON
    quitframe = sty.mySubframe(master)
    help = sty.myButton(quitframe, 'HELP', self.help_fun)
    help.grid(row=1, column=1)
    run_butt = sty.myButton(quitframe, 'RUN', self.nbrun)
    run_butt.grid(row = 1, column = 2)
    quit = sty.myButton(quitframe,'QUIT',self.quitall,bgc=tksty.qtcol)
    quit.grid(row = 1, column = 3)


    
  def init_dict(self):
    self.indict['shot']  = int(self.eshot.get().strip())
    self.indict['dt']    = float(self.edt  .get().strip())
    self.indict['tbeg']  = float(self.etbeg.get().strip())
    self.indict['tend']  = float(self.etend.get().strip())
    self.indict['n_rho'] = int(self.enrho.get().strip())
    self.indict['n_the'] = int(self.enthe.get().strip())
   
  def quitall(self):
#    sys.exit()
    mkroot.destroy()

  def error_dict(self):
    message_errordict = "Define DICT pressing the button 'dict' in the main window\n"
    tkMessageBox.showinfo("DICT error", message_errordict)

  def help_fun(self):
    message_general = "\n  GUI to generate input for nubeam_driver \n"
    message_info = "\n The data created will be written in ~/tr_client/TCV/<<SHOT NUMBER>>/\n"
    message_dict = "\n (1) Insert data needed in the frame on top and then press the button 'DICT'. This step is necessary otherwise the other functions won't work. \n"
    message_nml = "\n (2) Create the namelist using the button 'Namelist'. Before pressing 'Write Namelist', it is necessary to get the direction of current and magnetic field \n"
    message_1D  = "\n (3) Read and store 1D signals (Ip, Zeff, Vloop, NBI power)\n"
    message_2D  = "\n (4) Read and store 2D signals (ne, te, ti, vtor, p)\n"
    message_geo = "\n (5) Read and stores magnetic equilibrium (q, F, RZ)\n"
    message_run = "\n (6) Run NUBEAM with the data created now\n"
    message_authordate = "\n M. Vallar - 12/2016 (matteo.vallar@igi.cnr.it)\n"

    message = message_general+message_info+message_dict+message_nml+message_1D+message_2D+message_geo+message_run+message_authordate

    tkMessageBox.showinfo("HELPER", message)
    
  def geofread(self):
    try:
      self.indict['shot']
    except:
      self.error_dict()

    geo.GEO(self.indict)

  def dw1read(self):
    try:
      self.indict['shot']
      read_sig.DWR1(self.indict)
    except:
      self.error_dict()

  def dw2read(self):
    try:
      self.indict['shot']
      read_sig.DWR2(self.indict)
    except:
      self.error_dict()

  def write_nml(self):
    try:
      self.indict['shot']
    except:
      self.error_dict()

    nb_wrnml.WRNML(self.indict)

    
  def nbrun(self):
    try:
      self.indict['shot']
      nbrun.NBRUN(self.indict)
    except:
      #self.error_dict()
      nbrun.NBRUN(self.indict)


mkroot = tk.Tk(className=' NUBEAM input')
mkr = mkrun(mkroot)
mkroot.mainloop()
