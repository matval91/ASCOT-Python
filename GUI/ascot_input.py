#GUI
import Tkinter as tk
import tkFileDialog
# for the definition of appropriate plotting style
import tksty
import os

import ascot_prof, ascot_Bfield

sty = tksty.TKSTY()

def __init__():
    if __name__ == '__main__':
        fmframe = tk.TK()
    else:
        fmframe = tk.Toplevel()

    buttframe=sty.mySubframe(fmframe)
    sty.myButton(buttframe, 'B field', create_Binput, nrow=1, ncol=1)
    sty.myButton(buttframe, 'Profiles', create_profinput, nrow=1, ncol=2)
    #sty.myButton(buttframe, 'Input options', create_optinput, nrow=1, ncol=3)
    
    qframe = sty.mySubframe(fmframe)
    sty.myButton(qframe, 'QUIT', fmframe.destroy, bgc=tksty.qtcol)


def create_Binput():
#    if __name__ == '__main__':
#       fmframe = tk.TK()
#    else:
#        fmframe = tk.Toplevel()

    fname  = tkFileDialog.Open(filetypes=[('EQDSK','*.geq'),('All','*')],\
                               initialdir=os.environ['HOME'], title='Open file').show()
    if fname:      
        Bin = ascot_Bfield.Bfield_in(fname)

def create_profinput():
    if __name__ == '__main__':
       fmframe = tk.TK()
    else:
       fmframe = tk.Toplevel()

    #ASKS FOR nrho
    ion_dict={'Z':1, 'A':1, 'coll_mode':1}
    sty.myLabel(fmframe, 'Choose directory where to find input files').grid(row=1,column=1)
    dir_tmp = tkFileDialog.askdirectory(initialdir='/home/vallar', title='Choose dir with input files')
    sty.myLabel(fmframe,'# rho grid').grid(row=1, column=1)
    nrho_w = sty.myEntry(fmframe,119).grid(row=2, column=2)
    nrho = nrho_w.get()
    indict_frame = sty.mySubframe(fmframe)
    sty.myLabel(indict_frame,'Z').grid(row=1, column=1)
    sty.myLabel(indict_frame,'A').grid(row=1, column=3)
    sty.myLabel(indict_frame,'Collision mode').grid(row=2, column=1)

    Z_w = sty.myEntry(indict_frame,'1')
    A_w = sty.myEntry(indict_frame,'2')
    collmode_w= sty.myEntry(indict_frame,'1')
    
    Z_w.grid(row = 1, column = 2)
    A_w.grid    (row = 1, column = 4)
    collmode_w.grid   (row = 2, column = 2)
    ion_dict['Z'] = str(Z_w.get().strip())
    ion_dict['A'] = str(A_w.get().strip())
    ion_dict['coll_mode'] = str(collmode_w.get().strip())
    
    p = ascot_prof.profiles(dir_tmp, nrho, ion_dict)
    p.write_profs()
    
#def create_optinput():
#    if __name__ == '__main__':
#       fmframe = tk.TK()
#    else:
#       fmframe = tk.Toplevel()
       
    
    