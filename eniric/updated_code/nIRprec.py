"""
Created on Fri Feb  6 15:42:03 2015

@author: pfigueira
"""

import numpy as np
#import string

import matplotlib.pyplot as plt
from matplotlib import rc
#set stuff for latex usage
rc('text', usetex=True)

import IOmodule

import Qcalculator

atmmodel = "atmmodel/Average_TAPAS_2014.txt"
resampled_dir = "resampled_new/"
resampled_dir_OLD = "resampled/"


spectral_types = ["M0", "M3", "M6", "M9"]
bands = ["Y", "J", "H", "K"]
vsini = ["1.0", "5.0", "10.0"]
R = ["60k", "80k", "100k"]
sampling = ["3"]

def calculate_prec(plot_atm = False, plot_ste = False, paper_plots=True):
    
    print "Reading atmospheric model..."
    wav_atm, flux_atm, std_flux_atm, mask_atm = IOmodule.read_4col(atmmodel)

    wav_atm = np.array(wav_atm)/1000.0  #conversion from nanometers to micrometers
    flux_atm = np.array(flux_atm)
    std_flux_atm = np.array(std_flux_atm)
    mask_atm = np.array(mask_atm, dtype=bool)
    print("There were %d unmasked pixels out of %d." % (len(mask_atm[np.where(mask_atm)]), len(mask_atm)))
    print("The model range from %4.2f to %4.2f micron." % (wav_atm[0], wav_atm[-1]))
    print("Done.")
    
    if(plot_atm):
        
        #identify non-masked pixels
        selected_transmission = wav_atm[np.where(mask_atm)]
        
        plt.figure(1)
        plt.xlabel(r"wavelength [ $\mu$m ])")
        plt.ylabel(r"Transmission [ ] ")
        plt.plot(wav_atm, flux_atm, color ='k')
        plt.vlines(selected_transmission, 0.8, 1.0, colors = "b")
        plt.xlim(wav_atm[0], wav_atm[-1])
        plt.ylim(0.0, 1.0)
        #plt.legend(loc='best')
        plt.show()
        plt.close()        
    
    results = {}    #creating empty dictionary for the results
    
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        #print "Working on "+file_to_read+"."
                        wav_stellar, flux_stellar = IOmodule.read_2col(resampled_dir+file_to_read)
                        wav_stellar = np.array(wav_stellar)
                        flux_stellar = np.array(flux_stellar) / ((1.634e4)**2.0)       
                        
                        id_string = star+"-"+band+"-"+vel+"-"+resolution    # sample was left aside because only one value existed
                        
                        #gettin the wav, flux and mask from the atm model closes to the stellar

                        index_atm = np.searchsorted(wav_atm, wav_stellar)
                        #replace indexes outside the array, at the very end, by the value at the very end
                        index_atm = [index if(index<len(wav_atm)) else len(wav_atm)-1 for index in index_atm]   
                        
                        wav_atm_selected = wav_atm[index_atm]      
                        flux_atm_selected = flux_atm[index_atm]
                        mask_atm_selected = mask_atm[index_atm]                        
                        
                        if(id_string == "M0-J-1.0-100k"):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    #searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print "\tSanity Check: The S/N for the reference %s model was of %4.2f." % (id_string, SN_estimate) 
                            
                        if(plot_ste or plot_ste == id_string ):
                            #Plot the stellar spectr as considered 
                        
                            selected_transmission_stellar = wav_atm_selected[np.where(mask_atm_selected)]
                            
                            plt.figure(1)
                            plt.xlabel(r"wavelength [ $\mu$m ])")
                            plt.ylabel(r"Flux_stellar [ ] ")
                            plt.plot(wav_stellar, flux_stellar, color ='k')
                            plt.vlines(selected_transmission_stellar, np.min(flux_stellar), 0.3*np.max(flux_stellar), colors = "b")
                            plt.xlim(wav_stellar[0], wav_stellar[-1])
                            plt.ylim(np.min(flux_stellar)-0.1*(np.max(flux_stellar)-np.min(flux_stellar)), np.max(flux_stellar)+0.1*(np.max(flux_stellar)-np.min(flux_stellar)))
                            #plt.legend(loc='best')
                            plt.show()
                            plt.close()        
                            
                        #precision given by the first method:
                        prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)
                        
                        #precision as given by the second_method
                        """
                        Example Joao
                        a = np.array([1, 5, 6, 8, 16, 34, 5, 7, 10, 83, 12, 6, 17, 18])
                        b = np.array([1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1], dtype=bool)
                        
                        # this will give you a list of numpy arrays
                        c = np.array_split(a, np.where(np.diff(b))[0]+1)[::2]
                        
                        # this will give you a list of lists
                        d = [list(cc) for cc in c]
                        print d
                        >>> [[1, 5], [16, 34, 5], [83, 12], [17, 18]] 
                        """                        

                        wav_stellar_chunks_unformated = np.array_split(wav_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        wav_stellar_chunks = [list(chunk) for chunk in wav_stellar_chunks_unformated ]
                        
                        """
                        #test section
                        print "check that lengths are the same", len(wav_stellar), len(mask_atm_selected)
                        print("size of spectra %d vs number of chunks %d" % (len(wav_stellar), len(wav_stellar_chunks)))
                        print("number of true elements in all chunks: %d" % (len(mask_atm_selected[np.where(mask_atm_selected)])))
                        """
                        
                        flux_stellar_chunks_unformated = np.array_split(flux_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        flux_stellar_chunks = [list(chunk) for chunk in flux_stellar_chunks_unformated ]
                        
                        prec_2 = Qcalculator.RVprec_calc_chunks(wav_stellar_chunks, flux_stellar_chunks)                
                        
                        #precision as given by the third_method
                        prec_3 = Qcalculator.RV_prec_calc_Trans(wav_stellar, flux_stellar, flux_atm_selected)                  
                        
                        #adding result to the dictionary
                        results[id_string] = [prec_1, prec_2, prec_3]
    if(paper_plots):
        #print the paper plots
        """
        In every plot we will assume sample=3 and vsini=1
        y = RVprec between prec3 and prec2 
        x = different bands
        
        different panels for different spectral types
    
        different colors will represent the different resolutions
        """
        
        print("Results for vsini of 1.0 km/s")
        #preparation of data: plot1
        y1_60k = [results["M0-Y-1.0-60k"],results["M0-J-1.0-60k"], results["M0-H-1.0-60k"], results["M0-K-1.0-60k"]]        
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]        
        
        y1_80k = [results["M0-Y-1.0-80k"],results["M0-J-1.0-80k"], results["M0-H-1.0-80k"], results["M0-K-1.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]   
        
        y1_100k = [results["M0-Y-1.0-100k"],results["M0-J-1.0-100k"], results["M0-H-1.0-100k"], results["M0-K-1.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]        
        
        #preparation of data: plot2
        y2_60k = [results["M3-Y-1.0-60k"],results["M3-J-1.0-60k"], results["M3-H-1.0-60k"], results["M3-K-1.0-60k"]]        
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]        
        
        y2_80k = [results["M3-Y-1.0-80k"],results["M3-J-1.0-80k"], results["M3-H-1.0-80k"], results["M3-K-1.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]   
        
        y2_100k = [results["M3-Y-1.0-100k"],results["M3-J-1.0-100k"], results["M3-H-1.0-100k"], results["M3-K-1.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k] 
        
        #preparation of data: plot3
        y3_60k = [results["M6-Y-1.0-60k"],results["M6-J-1.0-60k"], results["M6-H-1.0-60k"], results["M6-K-1.0-60k"]]        
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]        
        
        y3_80k = [results["M6-Y-1.0-80k"],results["M6-J-1.0-80k"], results["M6-H-1.0-80k"], results["M6-K-1.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]   
        
        y3_100k = [results["M6-Y-1.0-100k"],results["M6-J-1.0-100k"], results["M6-H-1.0-100k"], results["M6-K-1.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k] 
        
        #preparation of data: plot4
        y4_60k = [results["M9-Y-1.0-60k"],results["M9-J-1.0-60k"], results["M9-H-1.0-60k"], results["M9-K-1.0-60k"]]        
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]        
        
        y4_80k = [results["M9-Y-1.0-80k"],results["M9-J-1.0-80k"], results["M9-H-1.0-80k"], results["M9-K-1.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]   
        
        y4_100k = [results["M9-Y-1.0-100k"],results["M9-J-1.0-100k"], results["M9-H-1.0-100k"], results["M9-K-1.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]         
        
        xticks = ["Y", "J", "H", "K"]        
        
        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])
                                   
        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)
        
        ax1.fill_between(range(1,5), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1,5), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1,5), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)
        
        ax1.text(1.0, positiony_max, "M0")        
        
        #ticks and labels
        ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(xticks)+0.5)    
        ax1.set_xticks(range(1,5))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))
        
        ax2 = fig.add_subplot(222)
        
        ax2.fill_between(range(1,5), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1,5), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1,5), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.text(1.0, positiony_max, "M3")           
        
        #ticks and labels
        ax2.set_xlim(0.5, len(xticks)+0.5)    
        ax2.set_xticks(range(1,5))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1,5), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1,5), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1,5), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.text(1.0, positiony_max, "M6")

        #ticks and labels
        ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')        
        ax3.set_xlim(0.5, len(xticks)+0.5)    
        ax3.set_xticks(range(1,5))
        ax3.set_xticklabels(xticks)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4 = fig.add_subplot(224)
        
        ax4.fill_between(range(1,5), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1,5), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1,5), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.text(1.0, positiony_max, "M9")
        
        #ticks and labels
        ax4.set_xlabel('Bands')        
        ax4.set_xlim(0.5, len(xticks)+0.5)    
        ax4.set_xticks(range(1,5))
        ax4.set_xticklabels(xticks)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        plt.show()
        plt.close()
        
        """
        VSINI=5km/s
        """
        print("Results for vsini of 5.0 km/s")
        
        #preparation of data: plot1
        y1_60k = [results["M0-Y-5.0-60k"],results["M0-J-5.0-60k"], results["M0-H-5.0-60k"], results["M0-K-5.0-60k"]]        
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]        
        
        y1_80k = [results["M0-Y-5.0-80k"],results["M0-J-5.0-80k"], results["M0-H-5.0-80k"], results["M0-K-5.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]   
        
        y1_100k = [results["M0-Y-5.0-100k"],results["M0-J-5.0-100k"], results["M0-H-5.0-100k"], results["M0-K-5.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]        
        
        #preparation of data: plot2
        y2_60k = [results["M3-Y-5.0-60k"],results["M3-J-5.0-60k"], results["M3-H-5.0-60k"], results["M3-K-5.0-60k"]]        
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]        
        
        y2_80k = [results["M3-Y-5.0-80k"],results["M3-J-5.0-80k"], results["M3-H-5.0-80k"], results["M3-K-5.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]   
        
        y2_100k = [results["M3-Y-5.0-100k"],results["M3-J-5.0-100k"], results["M3-H-5.0-100k"], results["M3-K-5.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k] 
        
        #preparation of data: plot3
        y3_60k = [results["M6-Y-5.0-60k"],results["M6-J-5.0-60k"], results["M6-H-5.0-60k"], results["M6-K-5.0-60k"]]        
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]        
        
        y3_80k = [results["M6-Y-5.0-80k"],results["M6-J-5.0-80k"], results["M6-H-5.0-80k"], results["M6-K-5.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]   
        
        y3_100k = [results["M6-Y-5.0-100k"],results["M6-J-5.0-100k"], results["M6-H-5.0-100k"], results["M6-K-5.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k] 
        
        #preparation of data: plot4
        y4_60k = [results["M9-Y-5.0-60k"],results["M9-J-5.0-60k"], results["M9-H-5.0-60k"], results["M9-K-5.0-60k"]]        
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]        
        
        y4_80k = [results["M9-Y-5.0-80k"],results["M9-J-5.0-80k"], results["M9-H-5.0-80k"], results["M9-K-5.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]   
        
        y4_100k = [results["M9-Y-5.0-100k"],results["M9-J-5.0-100k"], results["M9-H-5.0-100k"], results["M9-K-5.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]         
        
        xticks = ["Y", "J", "H", "K"]        
        
        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])
                                   
        fig = plt.figure(2)
        ax1 = fig.add_subplot(221)
        
        ax1.fill_between(range(1,5), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1,5), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1,5), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)
        
        ax1.text(1.0, positiony_max, "M0")        
        
        #ticks and labels
        ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(xticks)+0.5)    
        ax1.set_xticks(range(1,5))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))
        
        ax2 = fig.add_subplot(222)
        
        ax2.fill_between(range(1,5), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1,5), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1,5), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.text(1.0, positiony_max, "M3")           
        
        #ticks and labels
        ax2.set_xlim(0.5, len(xticks)+0.5)    
        ax2.set_xticks(range(1,5))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1,5), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1,5), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1,5), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.text(1.0, positiony_max, "M6")

        #ticks and labels
        ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')        
        ax3.set_xlim(0.5, len(xticks)+0.5)    
        ax3.set_xticks(range(1,5))
        ax3.set_xticklabels(xticks)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4 = fig.add_subplot(224)
        
        ax4.fill_between(range(1,5), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1,5), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1,5), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.text(1.0, positiony_max, "M9")
        
        #ticks and labels
        ax4.set_xlabel('Bands')        
        ax4.set_xlim(0.5, len(xticks)+0.5)    
        ax4.set_xticks(range(1,5))
        ax4.set_xticklabels(xticks)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        plt.show()
        plt.close()

        """
        VSINI=10km/s
        """
        print("Results for vsini of 10.0 km/s")

        #preparation of data: plot1
        y1_60k = [results["M0-Y-10.0-60k"],results["M0-J-10.0-60k"], results["M0-H-10.0-60k"], results["M0-K-10.0-60k"]]        
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]        
        
        y1_80k = [results["M0-Y-10.0-80k"],results["M0-J-10.0-80k"], results["M0-H-10.0-80k"], results["M0-K-10.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]   
        
        y1_100k = [results["M0-Y-10.0-100k"],results["M0-J-10.0-100k"], results["M0-H-10.0-100k"], results["M0-K-10.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]        
        
        #preparation of data: plot2
        y2_60k = [results["M3-Y-10.0-60k"],results["M3-J-10.0-60k"], results["M3-H-10.0-60k"], results["M3-K-10.0-60k"]]        
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]        
        
        y2_80k = [results["M3-Y-10.0-80k"],results["M3-J-10.0-80k"], results["M3-H-10.0-80k"], results["M3-K-10.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]   
        
        y2_100k = [results["M3-Y-10.0-100k"],results["M3-J-10.0-100k"], results["M3-H-10.0-100k"], results["M3-K-10.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k] 
        
        #preparation of data: plot3
        y3_60k = [results["M6-Y-10.0-60k"],results["M6-J-10.0-60k"], results["M6-H-10.0-60k"], results["M6-K-10.0-60k"]]        
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]        
        
        y3_80k = [results["M6-Y-10.0-80k"],results["M6-J-10.0-80k"], results["M6-H-10.0-80k"], results["M6-K-10.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]   
        
        y3_100k = [results["M6-Y-10.0-100k"],results["M6-J-10.0-100k"], results["M6-H-10.0-100k"], results["M6-K-10.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k] 
        
        #preparation of data: plot4
        y4_60k = [results["M9-Y-10.0-60k"],results["M9-J-10.0-60k"], results["M9-H-10.0-60k"], results["M9-K-10.0-60k"]]        
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]        
        
        y4_80k = [results["M9-Y-10.0-80k"],results["M9-J-10.0-80k"], results["M9-H-10.0-80k"], results["M9-K-10.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]   
        
        y4_100k = [results["M9-Y-10.0-100k"],results["M9-J-10.0-100k"], results["M9-H-10.0-100k"], results["M9-K-10.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]         
        
        xticks = ["Y", "J", "H", "K"]        
        
        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])
                                   
        fig = plt.figure(3)
        ax1 = fig.add_subplot(221)
        
        ax1.fill_between(range(1,5), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1,5), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1,5), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)
        
        ax1.text(1.0, positiony_max, "M0")        
        
        #ticks and labels
        ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(xticks)+0.5)    
        ax1.set_xticks(range(1,5))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))
        
        ax2 = fig.add_subplot(222)
        
        ax2.fill_between(range(1,5), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1,5), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1,5), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.text(1.0, positiony_max, "M3")           
        
        #ticks and labels
        ax2.set_xlim(0.5, len(xticks)+0.5)    
        ax2.set_xticks(range(1,5))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1,5), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1,5), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1,5), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.text(1.0, positiony_max, "M6")

        #ticks and labels
        ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')        
        ax3.set_xlim(0.5, len(xticks)+0.5)    
        ax3.set_xticks(range(1,5))
        ax3.set_xticklabels(xticks)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4 = fig.add_subplot(224)
        
        ax4.fill_between(range(1,5), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1,5), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1,5), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.text(1.0, positiony_max, "M9")
        
        #ticks and labels
        ax4.set_xlabel('Bands')        
        ax4.set_xlim(0.5, len(xticks)+0.5)    
        ax4.set_xticks(range(1,5))
        ax4.set_xticklabels(xticks)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        plt.show()
        plt.close()
        
    else:                    
        return results

###############################################################################    

def compare_runs():
    """
    Function that compares spectra as resampled in the two versions of the code
    """                   
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        #print "Working on "+file_to_read+"."
                        wav_stellar, flux_stellar = IOmodule.read_2col(resampled_dir+file_to_read)
                        wav_stellar = np.array(wav_stellar)
                        flux_stellar = np.array(flux_stellar) / ((1.634e4)**2.0)       
                        
                        wav_stellar_OLD, flux_stellar_OLD = IOmodule.read_2col(resampled_dir_OLD+file_to_read)
                        flux_stellar_OLD = np.array(flux_stellar_OLD) / ((1.634e4)**2.0)       

                        plt.figure(1)
                        plt.xlabel(r"wavelength [ $\mu$m ])")
                        plt.ylabel(r"Flux_stellar difference [ ] ")
                        plt.plot(wav_stellar, flux_stellar-flux_stellar_OLD, color ='k')
                        plt.show()
                        plt.close()     
                        
                        