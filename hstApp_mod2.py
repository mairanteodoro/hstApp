#!/usr/bin/env python

# NICE TUTORIAL: http://www.dreamincode.net/forums/topic/371440-tkinter-overview-with-a-fixed-width-grid/
import sys
import FileDialog

# resolving problem: Tkinter is for python 2; tkinter is for python 3
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkMessageBox, tkFileDialog
    import ttk
else:
    import tkinter as tk
    from tkinter import messagebox as tkMessageBox
    from tkinter import filedialog as tkFileDialog
    import tkinter.ttk as ttk

import numpy as np
import hstData, hstMeasure
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


class MainApplication(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        # defining the parent
        self.parent = parent
        self.parent.title('HST App - iso-velocity images module')

        # defining some containers
        self.mainContainer = tk.Frame(self.parent)
        self.mainContainer.pack(fill=tk.X, expand=1)
        tk.Grid.columnconfigure(self.mainContainer, 0, weight=1) # <- allows the button to expand to fill frame
        tk.Grid.columnconfigure(self.mainContainer, 1, weight=1) # <- allows the button to expand to fill frame
        tk.Grid.columnconfigure(self.mainContainer, 2, weight=1) # <- allows the button to expand to fill frame

        # controllers
        self.frame1 = tk.Frame(self.mainContainer, bg='red', bd=3, relief=tk.SUNKEN)
        self.frame1.grid(row=0, column=0, sticky=tk.W+tk.E+tk.N)
        tk.Grid.rowconfigure(self.frame1, 0, weight=1) # <- allows the button to expand to fill frame
        tk.Grid.columnconfigure(self.frame1, 0, weight=1) # <- allows the button to expand to fill frame
        tk.Grid.columnconfigure(self.frame1, 1, weight=1) # <- allows the button to expand to fill frame

        # time format
        self.frame2 = tk.Frame(self.mainContainer, bg='green', bd=3, relief=tk.SUNKEN)
        self.frame2.grid(row=1, column=0, sticky=tk.W+tk.E+tk.S)
        tk.Grid.columnconfigure(self.frame2, 0, weight=1)
        tk.Grid.columnconfigure(self.frame2, 1, weight=1)
        tk.Grid.columnconfigure(self.frame2, 2, weight=1)

        # canvas
        self.frame3 = tk.Frame(self.mainContainer, bg='blue', bd=3, relief=tk.SUNKEN)
        self.frame3.grid(row=0, column=1, columnspan=2, sticky='nsew')
        tk.Grid.rowconfigure(self.frame3, 0, weight=1)
        tk.Grid.columnconfigure(self.frame3, 0, weight=1)

        # parameters
        # self.frame3 = tk.Frame(self.parent, bg='blue', bd=3, relief=tk.SUNKEN)
        # self.frame3.pack(side=tk.BOTTOM, fill=tk.X, expand=1)
        # tk.Grid.rowconfigure(self.frame3, 0, weight=1)
        # tk.Grid.columnconfigure(self.frame3, 0, weight=1)


        # plot area
        # 1 - create a figure with no data...
        self.fig = Figure(figsize=(4,4))
        self.ax = self.fig.add_subplot(111)
        # 2 - ... create a canvas...
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame3)
        plt.gcf().canvas.draw()
        self.canvas.get_tk_widget().pack()
        # 3 - ... and finally a toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame3)
        self.toolbar.update()
        self.canvas._tkcanvas.pack()

        # call the widgets
        # self.plot_images_Button()
        # self.show_image_Button()
        # self.plot_flux_Button()
        # self.plot_position_Button()
        self.statusButton()
        self.quitButton()
        self.readDataButton()
        self.clearDataButton()
        self.velScale()
        # self.slitWidthScale()
        # self.radiusScale()
        # self.paScale()
        self.listLineIDListBox()
        self.logY_checkButton()
        self.timeFormatRadioButtons()

        # gridding and packing
        self.listLineIDListBox.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)

        self.statusButton.grid(row=1, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.quitButton.grid(row=1, column=1, sticky=tk.N+tk.S+tk.E+tk.W)

        self.readDataButton.grid(row=2, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.clearData.grid(row=2, column=1, sticky=tk.N+tk.S+tk.E+tk.W)

        # self.plot_images_Button.grid(row=3, columnspan=2, sticky=tk.N+tk.S+tk.E+tk.W)

        self.velLabel.grid(row=4, columnspan=2, sticky=tk.W+tk.E)
        self.velScale.grid(row=5, columnspan=2, sticky=tk.W+tk.E)

        # self.slitWidthLabel.grid(row=6, columnspan=2, sticky=tk.W+tk.E)
        # self.slitWidthScale.grid(row=7, columnspan=2, sticky=tk.W+tk.E)
        #
        # self.radiusLabel.grid(row=8, columnspan=2, sticky=tk.W+tk.E)
        # self.radiusScale.grid(row=9, columnspan=2, sticky=tk.W+tk.E)
        #
        # self.paLabel.grid(row=10, columnspan=2, sticky=tk.W+tk.E)
        # self.paScale.grid(row=11, columnspan=2, sticky=tk.W+tk.E)

        # self.logY_checkButton.pack()

        # Radio buttons for time format (called by self.timeFormatRadioButtons)
        self.timeFormatRadioButton1.grid(row=0, column=0, sticky=tk.W)
        self.timeFormatRadioButton2.grid(row=0, column=1, sticky=tk.W)
        self.timeFormatRadioButton3.grid(row=0, column=2, sticky=tk.W)
        self.timeFormatRadioButton4.grid(row=0, column=3, sticky=tk.W)

        # Button configurations
        ttk.Style().configure('green.TButton', foreground='green')
        ttk.Style().configure('blue.TButton', foreground='blue')
        ttk.Style().configure('red.TButton', foreground='red')


    # radio button options
    def printTimeSelection(self, var):
        if (var.get() == 1):
            # print('Time in JD')
            self.deleteCheckboxes()
            self.createCheckboxes(self.JD)
        if (var.get() == 2):
            from astropy.time import Time
            # print('Time in DecYear')
            self.deleteCheckboxes()
            date = np.asarray(self.JD, dtype=np.str)
            for i, item in enumerate(self.JD):
                time = Time(self.JD[i], format='jd')
                date[i] = '{0:.2f}'.format(time.decimalyear)
            self.createCheckboxes(date)
        if (var.get() == 3):
            # print('Time in Phi')
            self.deleteCheckboxes()
            phi = np.asarray(self.Phi, dtype=np.str)
            for i, item in enumerate(self.Phi):
                phi[i] = '{0:.4f}'.format(item)
            self.createCheckboxes(phi)
        if (var.get() == 4):
            from astropy.time import Time
            # print('Time in Greg')
            self.deleteCheckboxes()
            date = np.asarray(self.JD, dtype=np.str)
            for i, item in enumerate(self.JD):
                time = Time(self.JD[i], format='jd', out_subfmt='date')
                date[i] = '{0}'.format(time.iso)
            self.createCheckboxes(date)


    # Radio buttons for time format selection
    def timeFormatRadioButtons(self):
        var = tk.IntVar()
        self.timeFormatRadioButton1 = tk.Radiobutton(self.frame2, text="JD", variable=var, value=1, command=lambda:self.printTimeSelection(var))
        self.timeFormatRadioButton2 = tk.Radiobutton(self.frame2, text="DecYear", variable=var, value=2, command=lambda:self.printTimeSelection(var))
        self.timeFormatRadioButton3 = tk.Radiobutton(self.frame2, text="Phi", variable=var, value=3, command=lambda:self.printTimeSelection(var))
        self.timeFormatRadioButton4 = tk.Radiobutton(self.frame2, text="Greg", variable=var, value=4, command=lambda:self.printTimeSelection(var))


    ### Line ID listbox
    def listLineIDListBox(self):
        lineList = ('',) # <- create line (will contain self.lineID)
        self.listLineIDListBox = tk.Listbox(self.frame1, listvariable=tk.StringVar(value=lineList), width=5, height=5) # <- create list box

    def getLineIDFromListBox(self):
        # if there currently is nothing selected then return the last selected element
        if (len(self.listLineIDListBox.curselection()) == 0):
            return(self.listLineIDListBox.get(self.lastIndex))
        else:
            # get selected line index
            index = self.listLineIDListBox.curselection()[0]
            self.lastIndex = index
            # get the line's text
            return(self.listLineIDListBox.get(index))

    # show messages
    def printMessage(self):
        if (self.data != None):
            # if data exist enter here
            tkMessageBox.showinfo("HST App", "Data have been loaded from:\n\n" + self.fullPath)
        else:
            # warn user there is no data loaded
            tkMessageBox.showinfo("HST App", "\nData have not been loaded yet.")

    ### Check data status button
    def statusButton(self):
        self.statusButton = tk.Button(self.frame1, text='Status', command=self.printMessage)

    ### Quit button
    def quitButton(self):
        self.quitButton = tk.Button(self.frame1, text='Quit', command=self.confirmQuit)

    # confirm quitting
    def confirmQuit(self):
        answer = tkMessageBox.askyesno(title="HST App", message="\nDo you really want to quit?")
        if (answer):
            self.quit()

    ### Read data button
    def readDataButton(self):
        self.data = None # <- initialize data
        self.listEpochCheckButton = [] # <- initialize
        # self.readDataButton = tk.Button(self.frame1, text='Load data', command=self.readData)
        self.readDataButton = ttk.Button(self.frame1, text='Load data', command=self.readData, style='blue.TButton')

    # Return all checkboxes' status
    def printCheckButtons(self):
        for i, chk_bx in enumerate(self.listEpochCheckButton):
            self.listOfCheckboxes[i] = chk_bx.var.get()
        self.myCanvas() # <- call myCanvas() and plot the data associated to the selected checkboxes

    def deleteCheckboxes(self):
        for i, item in enumerate(self.listEpochCheckButton):
            self.listEpochCheckButton[i].grid_forget() # <- destroy checkboxes

    def createCheckboxes(self, time):
        # self.listOfCheckboxes = np.asarray(time) # number of checkboxes
        # creating the checkboxes with JD (epochs)
        self.listEpochCheckButton = []
        self.epoch_var = []
        col = 0
        row = 14
        for i, item in enumerate(time):
            var = tk.IntVar()
            self.epoch_checkButton = tk.Checkbutton(self.frame1,
                                                    relief=tk.RIDGE,
                                                    text=item,
                                                    variable=var,
                                                    onvalue=1, offvalue=0,
                                                    command=lambda:self.printCheckButtons(),
                                                    )
            self.epoch_checkButton.var = var
            self.listEpochCheckButton.append(self.epoch_checkButton)
            if (col >= 1):
                self.epoch_checkButton.grid(row=row, column=col, sticky=tk.W+tk.E)
                col = 0
                row += 1
            else:
                self.epoch_checkButton.grid(row=row, column=col, sticky=tk.W+tk.E)
                col += 1
            self.epoch_var.append(var)

    # reading data
    def readData(self):
        import os
        from os.path import expanduser
        home = expanduser("~")
        initialdir = home #'/Volumes/Kerberos/DATA/ETC/HST/TEDS_CUBE/NEW/'
        self.fullPath = dataList = tkFileDialog.askopenfilename(initialdir=initialdir,
                                                filetypes=[('JSON files', '*.json')])
        # dataDir = os.path.split(self.fullPath)[0]+'/'
        velMin = -60
        velMax = -20
        import hstData
        # for the App, dataList must have the full path to the FITS files
        self.data = hstData.read(dataList, velMin, velMax)
        # info about lines and epochs
        self.ntrans = len(np.unique(self.data['lineID'])) # number of lines
        self.nepoch = len(np.unique(self.data['JD'])) # number of observations
        self.JD = np.unique(self.data['JD']) # JD
        self.Phi = np.unique(self.data['Phi']) # orbital phase
        self.listOfCheckboxes = np.unique(self.data['JD']) # number of checkboxes
        # create checkboxes
        self.createCheckboxes(self.JD)
        # feed the list box with lineID
        lineList = np.unique(self.data['lineID'])
        self.listLineIDListBox.delete(0, tk.END)
        for i, item in enumerate(lineList):
            self.listLineIDListBox.insert(tk.END, item)
        self.lineID = self.listLineIDListBox.get(0)
        # warn the user the data have been loaded
        tkMessageBox.showinfo("HST App", "\nData have been succesfully loaded from '" + os.path.split(self.fullPath)[1] + "'")
        self.listLineIDListBox.select_set(0) #This only sets focus on the first item.
        self.listLineIDListBox.event_generate("<<ListboxSelect>>")
        self.lastIndex = None

    ### Clear data from current session
    def clearDataButton(self):
        self.clearData = ttk.Button(self.frame1, text='Clear data', command=self.confirmClearData, style='red.TButton')

    # confirm clearing data
    def confirmClearData(self):
        answer = tkMessageBox.askyesno(title="HST App", message="\nAre you sure you want to clear the loaded data?")
        if (answer):
            self.data = None # <- clear data
            self.listLineIDListBox.delete(0, tk.END) # <- clear list box
            # delete checkboxes
            self.deleteCheckboxes()
            # clear plot area
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            self.fig.canvas.draw()
            # warn the user everything is clear
            tkMessageBox.showwarning(title="HST App", message="\nData has been deleted.")

    ### Velocity scale
    def velScale(self):
        self.velVar = tk.StringVar()
        self.velLabel = tk.Label(self.frame1, textvariable=self.velVar)
        self.velScale = tk.Scale(self.frame1, from_=-500, to=+500,
                            orient=tk.HORIZONTAL, resolution=20,
                            sliderlength=20, showvalue=0,
                            length=200, width=20,
                            command=self.onVelScale)
    # update velLabel
    # note: 'val' is the new scale value (a str) passed in by tk.Scale() every time the slider moves
    def onVelScale(self, val):
        self.velVar.set("Velocity: {:+0.0f} km/s".format(float(val)))

    ### Slit width scale
    def slitWidthScale(self):
        self.slitWidthVar = tk.StringVar()
        self.slitWidthLabel = tk.Label(self.frame1, textvariable=self.slitWidthVar)
        self.slitWidthScale = tk.Scale(self.frame1, from_=0.1, to=10,
                            orient=tk.HORIZONTAL, resolution=0.1,
                            sliderlength=20, showvalue=0,
                            length=200, width=20,
                            command=self.onSlitWidthScale)
    def onSlitWidthScale(self, val):
        self.slitWidth = float(val)
        self.slitWidthVar.set("Slit width: {:0.1f} arcsec".format(float(val)))

    ### Radius scale
    def radiusScale(self):
        self.radiusVar = tk.StringVar()
        self.radiusLabel = tk.Label(self.frame1, text="Velocity:", textvariable=self.radiusVar)
        self.radiusScale = tk.Scale(self.frame1, from_=0.1, to=10,
                            orient=tk.HORIZONTAL, resolution=0.1,
                            sliderlength=20, showvalue=0,
                            length=200, width=20,
                            command=self.onRadiusScale)
    def onRadiusScale(self, val):
        self.radius = float(val)
        self.radiusVar.set("Radius: {:0.1f} arcsec".format(float(val)))

    ### PA scale
    def paScale(self):
        self.paVar = tk.StringVar()
        self.paLabel = tk.Label(self.frame1, text="Velocity:", textvariable=self.paVar)
        self.paScale = tk.Scale(self.frame1, from_=0, to=359,
                            orient=tk.HORIZONTAL, resolution=1,
                            sliderlength=20, showvalue=0,
                            length=200, width=20,
                            command=self.onPaScale)
    def onPaScale(self, val):
        self.pa = float(val)
        self.paVar.set(u"PA: {:0.1f}\N{DEGREE SIGN} ".format(float(val)))

    ### Change Y scale of plot
    def logY_checkButton(self):
        self.logY_checkButton_var = tk.StringVar()
        self.logY_checkButton = tk.Checkbutton(self.frame3,
                                                text="log10/linear",
                                                variable=self.logY_checkButton_var,
                                                onvalue="log", offvalue="linear",
                                                command=lambda:self.myCanvas())
        self.logY_checkButton.deselect()

    ### Plot data button -> canvas
    def plot_images_Button(self):
        self.counter = 0
        self.plot_images_Button = ttk.Button(self.frame1, text='Plot radial profile', command=self.myCanvas, style='green.TButton')

    # plot data
    def plot_images(self, x, y, epoch):
        if (self.data != None):
            # if data exist enter here
            self.fig.clf() # <- clear the entire figure instance
            self.ax = self.fig.add_subplot(111) # <- create new axes
            self.ax.set_yscale(self.logY_checkButton_var.get())
            self.ax.set_title(self.lineID)
            for i, item in enumerate(x):
                self.ax.plot(x[i], y[i],
                            markersize=3,
                            linewidth=1,
                            label='{0:.3f}'.format(epoch[i]),
                            ) # <- plot on new axes
            self.ax.legend(loc='best', ncol=2, fontsize=8) # <- add legend
            self.fig.canvas.draw() # <- draw figure with new axes
        else:
            # warn user there is no data loaded
            self.printMessage()

    def isoVelImage(self, data, velocity=None):
            ionData_prof1 = [] # np.empty((len(lineID), nepoch))
            ionData_x1 = [] # np.empty((len(lineID), nepoch))
            profile = hstMeasure.profile(data, radius, pa, slitWidth)
            ionData1 = hstData.get(profile, key='lineID', value=self.lineID) # retrieving the data for the ion
            ionData_prof1.extend([[ionData1[x]['profile'],ionData1[x]['lineID']] for x in range(self.nepoch)]) # storing the distance
            ionData_x1.extend([[ionData1[x]['radPos'],ionData1[x]['lineID']] for x in range(self.nepoch)]) # storing PA
            return([ionData_x1, ionData_prof1])

    ### Canvas
    def myCanvas(self):
        import numpy as np
        if (self.data != None):
            # if data exist enter here
            self.lineID = self.getLineIDFromListBox()
            initial_slitWidth = float(self.slitWidth)
            initial_radius = float(self.radius)
            initial_pa = float(self.pa)
            res = self.isoVelImage(self.data,
                                    velocity=initial_slitWidth,
                                    )
            # selected epochs (checkboxes = 1)
            nonZero = np.flatnonzero(self.listOfCheckboxes)
            x1, prof1 = res[0], res[1]
            xvar, yvar = [], []
            for i, item in enumerate(nonZero):
                xvar.append(x1[item][0])
                yvar.append(prof1[item][0])
            # img = ionData[0]['image']
            # print("Hello, from inside 'canvas'!")
            epoch = self.JD[nonZero]
            self.plot_images(xvar, yvar, epoch)
        else:
            # there is no data loaded
            a, b, size = 0., 360., 101
            xvar = np.linspace(a, b, num=size, endpoint=True)
            yvar = np.sin(xvar / 180. * np.pi) #(b - a) * np.random.random_sample(size) + a
            self.plot_images(xvar, yvar)



if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1200x700+10+10") # format: 'wxh+/-x+/-y'
    root.resizable(0, 0) # disabling window resizing
    MainApplication(root).pack()
    root.mainloop()


#
# from http://stackoverflow.com/questions/17466561/best-way-to-structure-a-tkinter-application/17470842#17470842
#
# class Navbar(tk.Frame): ...
# class Toolbar(tk.Frame): ...
# class Statusbar(tk.Frame): ...
# class Main(tk.Frame): ...
#
# class MainApplication(tk.Frame):
#     def __init__(self, parent, *args, **kwargs):
#         tk.Frame.__init__(self, parent, *args, **kwargs)
#         self.statusbar = Statusbar(self, ...)
#         self.toolbar = Toolbar(self, ...)
#         self.navbar = Navbar(self, ...)
#         self.main = Main(self, ...)
#
#         self.statusbar.pack(side="bottom", fill="x")
#         self.toolbar.pack(side="top", fill="x")
#         self.navbar.pack(side="left", fill="y")
#         self.main.pack(side="right", fill="both", expand=True)
