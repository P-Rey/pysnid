# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 09:53:14 2018

@author: Peter
"""

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
## Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import tkinter as tk
from matplotlib import style
import matplotlib.backends.backend_tkagg as tkagg

    
class SeaofBTCapp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
#        tk.Tk.iconbitmap(self, default ="@/Downloads/ucd_brandmark_colour.ico")
        tk.Tk.wm_title(self, "PySNID")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}



        frame = StartPage(container, self)

        self.frames[StartPage] = StartPage(container, self) # frame

        frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()


def qf(quickPrint):
    print(quickPrint)

        
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)

        root = tk.Tk()
    
       # root.wm_title("PySNID")
        fig = Figure(figsize=(10, 4), dpi=100)
        fig1 = Figure(figsize=(10, 4), dpi=100)

        fig.add_subplot(111).plot(input_spectrum[:,0],input_spectrum[:,1])
        
        fig.add_subplot(111).plot(template_spectrum[:,0],template_spectrum[:,1])
        fig1.add_subplot(121).plot(correlation)
        
        canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        tkagg.NavigationToolbar2Tk(canvas, root)
        
        canvas1 = FigureCanvasTkAgg(fig1, master=root)  # A tk.DrawingArea.
        canvas1.draw()
        canvas1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        
        # canvas is your canvas, and root is your parent (Frame, TopLevel, Tk instance etc.)
        tkagg.NavigationToolbar2Tk(canvas1, root)
        
        tk.Label(root, text=rlap%(0,0), borderwidth=1 ).grid(row=0,column=0)


if __name__ == "__main__":
    
    import reader
    template_spectrum = reader.readspec('SN2011fe_2011-08-25_00-00-00_TNG_DOLORES_PTF.ascii')
    input_spectrum = reader.parsereadspec()
    import apodise.apodiser as ap
    apodised_input_spectrum = ap.apodise(input_spectrum)
    apodised_template_spectrum = ap.apodise(template_spectrum)  
    import logtrim
    input_limits = logtrim.getspeclimits(apodised_input_spectrum)
    template_limits = logtrim.getspeclimits(apodised_template_spectrum)
    log_input_spectrum = logtrim.logspec(apodised_input_spectrum)
    log_template_spectrum = logtrim.logspec(apodised_template_spectrum)
    log_input_spectrum, log_template_spectrum = logtrim.logspectrim(log_input_spectrum, log_template_spectrum, input_limits, template_limits)
    import correlator
    correlation, rlap = correlator.correlate_spec(log_input_spectrum, log_template_spectrum)
    
    
    LARGE_FONT= ("Verdana", 12)
    app = SeaofBTCapp()
    app.mainloop()