from __future__ import division

import sys

if __name__ == '__main__':

    # This code computes k(lam) using the Calzetti curve.

    rv = 4.05

    # Define line wavelngths. In air.
    line_name_list = ['Hbeta', '[OIII]4959', '[OIII]5007', '[OI]6300', '[OI]6364', '[NII]6548', 'Halpha', '[NII]6583', '[SII]6716', '[SII]6731'] 
    line_wav_list  = [0.4861, 0.4959, 0.5007, 0.63, 0.6364, 0.6548, 0.6563, 0.65835, 0.6716, 0.6731]  # has to be in microns

    for i in range(len(line_wav_list)):

        line_wav = line_wav_list[i]

        if line_wav > 0.63:
            klam = 2.659 * (-1.857 + 1.040/line_wav) + rv

        elif line_wav < 0.63:
            klam = 2.659 * (-2.156 + 1.509/line_wav - 0.198/line_wav**2 + 0.011/line_wav**3) + rv

        elif line_wav == 0.63:
            klam = 3.49
            # Since the curves dont exactly meet at 0.63 micron, I've taken the average of the two.
            # They're close though. One gives me klam=3.5 and the other klam=3.48

        print line_name_list[i], "{:.4}".format(klam)