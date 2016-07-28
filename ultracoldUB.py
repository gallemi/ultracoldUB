import sys
import subprocess
from contextlib import contextmanager
import os

#c_bs = 0
#c_ds = 0

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

while True:

    print("\nChoose one case:")
    print("\t(1) Wavepacket dispersion")
    print("\t(2) Dark solitons")
    print("\t(3) Bright solitons")
    print("Or type 'q' to quit.")

    case = raw_input(" > ")

    try:
        int(case)

    except ValueError:
        try:
            str(case)
        except ValueError:
            continue
        else:
            if case == "q" or case == "Q":
                sys.exit("End of program")
            else:
                print("Not an available option.")
                continue

    else:

        if int(case) == 1:
            print("Wavepacket dispersion")
            with cd('./Wavepackdisper'):
                os.system('python gpe_fft_ts_WP_v1.py')

        elif int(case) == 2:
            print("Dark solitons")
            with cd('./darksolitons'):
                os.system('python gpe_fft_ts_DS_v1.py')
            #if c_ds == 0:
            #    import darksolitons.gpe_fft_ts_DS_v1 as ds
            #    c_ds = 1
            #else:
            #    reload(ds)

        elif int(case) == 3:
            print("Bright solitons")
            with cd('./brightsolitons'):
                os.system('python gpe_bright_solitons.py')
            #if c_bs == 0:
            #    import brightsolitons.gpe_bright_solitons as bs
            #    c_bs = 1
            #else:
            #    reload(bs)

        else:
            print("Number must be 1, 2 or 3.")
            continue
