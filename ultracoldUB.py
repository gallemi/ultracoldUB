import sys
import subprocess
from contextlib import contextmanager
import os
try:
    from tkinter import *  
    version='3'
except ImportError:
    from Tkinter import *
    version=''

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)



class Demomain:
    def __init__(self,master):
        self.name=""
        self.master=master
        frame=Frame(self.master)
        frame.pack()
        self.v=IntVar()
        options = [
            ("Wave packet dispersion",1),
            ("Dark Soliton",2),
            ("Bright Soliton",3)
        ]
        
        label = Label(frame, text="Choose one case:", background="black", foreground="white",font = "Verdana 16 bold")
        label.pack(fill=X)
        for txt, val in options:
            Radiobutton(master, 
                text=txt,
                variable=self.v, 
                command=self.show_values,
                value=val).pack(anchor=W)
        self.button=Button(frame, text='OK', command=self.exitroot).pack()


    def show_values(self):
        self.name= (self.v.get())
        
    def exitroot(self):
        root.destroy()
        
root = Tk()
root.wm_title("ULTRACOLD UB")
Demain=Demomain(root)
root.mainloop()

case=Demain.name


#while True:
#
#    print("\nChoose one case:")
#    print("\t(1) Wavepacket dispersion")
#    print("\t(2) Dark solitons")
#    print("\t(3) Bright solitons")
#    print("Or type 'q' to quit.")
#
#    case = raw_input(" > ")
#
#    try:
#        int(case)
#
#    except ValueError:
#        try:
#            str(case)
#        except ValueError:
#            continue
#        else:
#            if case == "q" or case == "Q":
#                sys.exit("End of program")
#            else:
#                print("Not an available option.")
#                continue
#
#    else:

if int(case) == 1:
    print("Wavepacket dispersion")
    with cd('./Wavepackdisper'):
        os.system('python'+version+' gpe_fft_ts_WP_v1.py')

elif int(case) == 2:
    print("Dark solitons")
    with cd('./darksolitons'):
        os.system('python'+version+' gpe_fft_ts_DS_v1.py')

elif int(case) == 3:
    print("Bright solitons")
    with cd('./brightsolitons'):
        os.system('python gpe_bright_solitons.py')
#
#        else:
#            print("Number must be 1, 2 or 3.")
#            continue
