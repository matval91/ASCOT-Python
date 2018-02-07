import sys
# GUI 
import Tkinter as tk
import tksty
import tkMessageBox

#ascot-related
import ascot_output
#other stuff
import os,fnmatch,shutil
#import matplotlib.pylab as plt

from os.path import expanduser

sty=tksty.TKSTY()

class mkrun:

  def __init__(self, master):
      # FRAME FOR NAME, author and date
      nameframe=sty.mySubframe(master)

      # FRAME FOR DICT WITH MACHINE,RUN, SHOT
      inputframe=sty.mySubframe(master)
      
      sty.myLabel(inputframe,'Machine').grid(row=1, column=1)
      sty.myLabel(inputframe,'Run'    ).grid(row=1, column=3)
      sty.myLabel(inputframe,'Shot'   ).grid(row=2, column=1)

      machine_w = sty.myEntry(inputframe,'JT-60SA')
      run_w     = sty.myEntry(inputframe,'002'    )
      shot_w    = sty.myEntry(inputframe,'111'    )
    
      machine_w.grid(row = 1, column = 2)
      run_w.grid    (row = 1, column = 4)
      shot_w.grid   (row = 2, column = 2)
      self.old_dir=os.getcwd()
      machine = str(machine_w.get().strip())
      run     = str(run_w.get().strip())
      shot    = str(shot_w.get().strip())
      home=os.path.expanduser("~")
      dir=home+'/ASCOT/runs/'+str(machine)+'/'+str(run)+'/'+str(shot)+'/'
      sty.myButton(inputframe,'DIR', os.chdir(dir), nrow=2, ncol=4)

      # FRAME for input/output
      inoutframe = sty.mySubframe(master)
      #in_but  = sty.myButton(inoutframe,'INPUT',self.input)
      #in_but.grid(row=1, column=1)
      out_but = sty.myButton(inoutframe,'OUTPUT', self.output) 
      out_but.grid(row=1, column=2)
      
      # FRAME FOR QUIT BUTTON
      quitframe = sty.mySubframe(master)
      help = sty.myButton(quitframe, 'HELP', self.help_fun)
      help.grid(row=1, column=1)
      quit = sty.myButton(quitframe,'QUIT',self.quitall,bgc=tksty.qtcol)
      quit.grid(row = 1, column = 2)


  #def input(self):

  def output(self):
    print os.getcwd()
    ascot_output.main()


  def quitall(self):
    os.chdir(self.old_dir)
    mkroot.destroy()

  def help_fun(self):
    message_general = "\n  GUI to look at output or generate input for ASCOT and BBNBI \n"
    message_info = "\n The data created will be written in ~/ASCOT/run/<<MACHINE>>/<<RUN>>/<<SHOT>> \n"
    message_dir = "\n (1) Insert data needed in the frame on top and then press the button 'DIR'. This step is necessary otherwise the other functions won't use the data you want. \n"
    message_authordate = "\n M. Vallar - 2/2017 (matteo.vallar@igi.cnr.it)\n"

    message = message_general+message_info+message_dir+message_authordate

    tkMessageBox.showinfo("HELPER", message)


mkroot = tk.Tk(className='Ascot GUI')
mkr = mkrun(mkroot)
mkroot.mainloop()
