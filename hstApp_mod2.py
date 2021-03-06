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

import pdb as pdb


class MainApplication(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        # defining the parent
        self.parent = parent
        self.parent.title('HST App - iso-velocity images module')

        # defining some containers
        self.mainContainer = tk.Frame(self.parent)
        self.mainContainer.pack(fill=tk.X, expand=1)
        # tk.Grid.columnconfigure(self.mainContainer, 0, weight=1) # <- allows the button to expand to fill frame
        # tk.Grid.columnconfigure(self.mainContainer, 1, weight=1) # <- allows the button to expand to fill frame
        # tk.Grid.columnconfigure(self.mainContainer, 2, weight=1) # <- allows the button to expand to fill frame

        # controllers
        # self.frame1 = tk.Frame(self.mainContainer, bg='red', bd=3, relief=tk.SUNKEN)
        self.frame1 = tk.Frame(self.mainContainer, bg='red', bd=3, relief=tk.SUNKEN)
        self.frame1.grid(row=0, column=0, sticky=tk.W+tk.E+tk.N)
        # tk.Grid.rowconfigure(self.frame1, 0, weight=1) # <- allows the button to expand to fill frame
        # tk.Grid.columnconfigure(self.frame1, 0, weight=1) # <- allows the button to expand to fill frame
        # tk.Grid.columnconfigure(self.frame1, 1, weight=1) # <- allows the button to expand to fill frame

        # time format
        # self.frame2 = tk.Frame(self.mainContainer, bg='green', bd=3, relief=tk.SUNKEN)
        self.frame2 = tk.Frame(self.mainContainer, bg='green', bd=3, relief=tk.SUNKEN)
        self.frame2.grid(row=1, column=0, sticky=tk.S)
        # tk.Grid.columnconfigure(self.frame2, 0, weight=1)
        # tk.Grid.columnconfigure(self.frame2, 1, weight=1)
        # tk.Grid.columnconfigure(self.frame2, 2, weight=1)
        # tk.Grid.columnconfigure(self.frame2, 3, weight=1)

        # canvas
        # self.frame3 = tk.Frame(self.mainContainer, bg='blue', bd=3, relief=tk.SUNKEN)
        self.frame3 = tk.Frame(self.mainContainer, bg='blue', bd=3, relief=tk.SUNKEN)
        self.frame3.grid(row=0, column=1, columnspan=3, sticky='nsew')
        # tk.Grid.rowconfigure(self.frame3, 0, weight=1)
        # tk.Grid.columnconfigure(self.frame3, 0, weight=1)

        # display options controllers
        # self.frame4 = tk.Frame(self.mainContainer, bg='blue', bd=3, relief=tk.SUNKEN)
        self.frame4 = tk.Frame(self.mainContainer, bg='purple', bd=3, relief=tk.SUNKEN)
        self.frame4.grid(row=1, column=1, sticky=tk.W)
        # tk.Grid.rowconfigure(self.frame4, 1, weight=1)
        # tk.Grid.columnconfigure(self.frame4, 0, weight=1)
        # tk.Grid.columnconfigure(self.frame4, 1, weight=1)
        # tk.Grid.columnconfigure(self.frame4, 2, weight=1)

        # re-scale option
        self.frame5 = tk.Frame(self.mainContainer, bg='orange', bd=3, relief=tk.SUNKEN)
        self.frame5.grid(row=1, column=1, sticky=tk.E)

        # parameters
        # self.frame3 = tk.Frame(self.parent, bg='blue', bd=3, relief=tk.SUNKEN)
        # self.frame3.pack(side=tk.BOTTOM, fill=tk.X, expand=1)
        # tk.Grid.rowconfigure(self.frame3, 0, weight=1)
        # tk.Grid.columnconfigure(self.frame3, 0, weight=1)


        # plotting area
        # 1 - create a figure with no data...
        self.fig = Figure(figsize=(4,4))
        self.ax = self.fig.add_subplot(111)
        # 2 - ... create a canvas...
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame3)
        plt.gcf().canvas.draw()
        self.ax.axes.set_visible(False)
        self.canvas.get_tk_widget().pack()
        # 3 - ... and finally a toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame3)
        self.toolbar.update()
        self.canvas._tkcanvas.pack()


        # call the widgets
        self.plot_images_Button()
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
        self.timeFormatRadioButtons()
        self.displayOptionsRadioButtons()
        self.rescale_checkButton()

        # gridding and packing
        self.listLineIDListBox.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)

        self.statusButton.grid(row=1, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.quitButton.grid(row=1, column=1, sticky=tk.N+tk.S+tk.E+tk.W)

        self.readDataButton.grid(row=2, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.clearData.grid(row=2, column=1, sticky=tk.N+tk.S+tk.E+tk.W)

        self.plot_images_Button.grid(row=3, columnspan=2, sticky=tk.N+tk.S+tk.E+tk.W)

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

        self.rescale_checkButton.pack()

        # Radio buttons for time format (called by self.timeFormatRadioButtons)
        self.timeFormatRadioButton1.grid(row=0, column=0, sticky=tk.W)
        self.timeFormatRadioButton2.grid(row=0, column=1, sticky=tk.W)
        self.timeFormatRadioButton3.grid(row=0, column=2, sticky=tk.W)
        self.timeFormatRadioButton4.grid(row=0, column=3, sticky=tk.W)

        # Radio buttons for display options
        self.displayOptionsRadioButton1.grid(row=0, column=0, sticky=tk.W)
        self.displayOptionsRadioButton2.grid(row=0, column=1, sticky=tk.W)
        self.displayOptionsRadioButton3.grid(row=0, column=2, sticky=tk.W)

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
        self.timeFormatVar = tk.IntVar()
        self.timeFormatVar.set(1)
        self.timeFormatRadioButton1 = tk.Radiobutton(self.frame2, text="JD", variable=self.timeFormatVar, value=1, command=lambda:self.printTimeSelection(self.timeFormatVar))
        self.timeFormatRadioButton2 = tk.Radiobutton(self.frame2, text="DecYear", variable=self.timeFormatVar, value=2, command=lambda:self.printTimeSelection(self.timeFormatVar))
        self.timeFormatRadioButton3 = tk.Radiobutton(self.frame2, text="Phi", variable=self.timeFormatVar, value=3, command=lambda:self.printTimeSelection(self.timeFormatVar))
        self.timeFormatRadioButton4 = tk.Radiobutton(self.frame2, text="Greg", variable=self.timeFormatVar, value=4, command=lambda:self.printTimeSelection(self.timeFormatVar))


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
        self.dataList = None
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
        # creating the checkboxes with time labels
        self.listEpochCheckButton = []
        self.epoch_var = time
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
        pdb.set_trace()
        # self.epoch_checkButton.select()
        # self.epoch_var.append(var)

    def myRescale(self):
        # pdb.set_trace()
        if (self.rescale_checkButton_var.get() == '0'):
            # no re-scaling
            img1 = self.isoVelImage()
            img = img1
        if (self.rescale_checkButton_var.get() == '1'):
            # re-scaling
            import scipy.ndimage
            img2 = self.isoVelImage()
            img_ref = img2[0][0] # <- reference image
            img = img2
            for i in range(len(img)):
                factor = (img[i][0]).shape[0] / (img_ref).shape[0]
                pdb.set_trace()
                img[i][0] = (img2[i][0])
        # call plotting function
        self.plot_images(img)

    ### Re-scale images to match each other
    def rescale_checkButton(self):
        self.rescale_checkButton_var = tk.StringVar()
        self.rescale_checkButton = tk.Checkbutton(self.frame5,
                                                text="Re-scale images",
                                                variable=self.rescale_checkButton_var,
                                                onvalue=1, offvalue=0,
                                                command=lambda:self.myRescale())
        self.rescale_checkButton.deselect()

    # reading data
    def readData(self):
        import os
        from os.path import expanduser
        home = expanduser("~")
        initialdir = home #'/Volumes/Kerberos/DATA/ETC/HST/TEDS_CUBE/NEW/'
        self.fullPath = self.dataList = tkFileDialog.askopenfilename(initialdir=initialdir,
                                                filetypes=[('JSON files', '*.json')])
        # dataDir = os.path.split(self.fullPath)[0]+'/'
        print("++++++")
        print(self.dataList)
        print("++++++")
        velMin = -40#float(self.curVel)
        velMax = -40#float(self.curVel)
        import hstData
        # for the App, dataList must have the full path to the FITS files
        self.data = hstData.read(self.dataList, velMin, velMax)
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
        tkMessageBox.showinfo("HST App", "\nData have been succesfully loaded from '" + os.path.split(self.dataList)[1] + "'")
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
        self.curVel = float(val)
        self.velVar.set("Velocity: {:+0.0f} km/s".format(float(val)))
        if (self.data != None):
            self.myCanvas()

    # display options functionality
    def displayOptions(self, var):
        # pdb.set_trace()
        if (var.get() == 1):
            # linear scale
            img1 = self.isoVelImage()
            img = img1
        if (var.get() == 2):
            # log10 scale
            img2 = self.isoVelImage()
            img = img2
            for i in range(len(img)):
                img[i][0] = np.log10(img2[i][0])
        if (var.get() == 3):
            # square root scale
            img3 = self.isoVelImage()
            img = img3
            for i in range(len(img)):
                img[i][0] = np.sqrt(img3[i][0])
        # call plotting function
        self.plot_images(img)

    # display options radio buttons
    def displayOptionsRadioButtons(self):
        self.displayOption = tk.IntVar()
        self.displayOption.set(1)
        self.displayOptionsRadioButton1 = tk.Radiobutton(self.frame4, text="linear", variable=self.displayOption, value=1, command=lambda:self.displayOptions(self.displayOption))
        self.displayOptionsRadioButton2 = tk.Radiobutton(self.frame4, text="log10", variable=self.displayOption, value=2, command=lambda:self.displayOptions(self.displayOption))
        self.displayOptionsRadioButton3 = tk.Radiobutton(self.frame4, text="square root", variable=self.displayOption, value=3, command=lambda:self.displayOptions(self.displayOption))

    ### Plot data button -> canvas
    def plot_images_Button(self):
        self.counter = 0
        self.plot_images_Button = ttk.Button(self.frame1, text='Update display', command=self.myCanvas, style='green.TButton')

    # plot data
    def plot_images(self, img):
        if (self.data != None):
            # if data exist enter here
            self.fig.clf() # <- clear the entire figure instance
            self.ax = self.fig.add_subplot(111) # <- create new axes
            self.ax.set_xlabel(r'$\Delta\alpha$ (arcsec)', fontsize=10)
            self.ax.set_ylabel(r'$\Delta\delta$ (arcsec)', fontsize=10)
            self.ax.tick_params(axis='both', which='major', labelsize=10)
            # pdb.set_trace()
            for i in range(len(self.nonZero)):
                self.ax.set_title('{0} ({1})'.format(self.lineID, self.selectedEpochs[i]))
                self.ax.imshow(img[self.nonZero[i]][0], # <- plot on new axes
                                interpolation='none',
                                extent=[img[self.nonZero[i]][2][0]-(-0.5*img[self.nonZero[i]][4]),
                                        img[self.nonZero[i]][2][-1]+(-0.5*img[self.nonZero[i]][4]),
                                        img[self.nonZero[i]][3][0]+(-0.5*img[self.nonZero[i]][4]),
                                        img[self.nonZero[i]][3][-1]-(-0.5*img[self.nonZero[i]][4])],
                                )
            self.fig.canvas.draw() # <- draw figure with new axes
        else:
            # warn user there is no data loaded
            self.printMessage()

    def isoVelImage(self):
            ionData_img = [] # np.empty((len(lineID), nepoch))
            import hstData
            # for the App, dataList must have the full path to the FITS files
            data = hstData.read(self.dataList, self.curVel, self.curVel)
            ionData = hstData.get(data, key='lineID', value=self.lineID) # retrieving the data for the ion
            ionData_img.extend([[ionData[x]['image'],
                                ionData[x]['lineID'],
                                ionData[x]['xScale'],
                                ionData[x]['yScale'],
                                ionData[x]['pixScale']] for x in range(self.nepoch)]) # storing image, lineID, xCenterOfPixels, yCenterOfPixels, & pixelScale
            return(ionData_img)

    ### Canvas
    def myCanvas(self):
        # import numpy as np
        if (self.data != None):
            # if data exist enter here
            self.lineID = self.getLineIDFromListBox()
            initial_vel = float(self.curVel)
            img = self.isoVelImage()
            # index of selected epochs (checkboxes = 1)
            self.nonZero = np.flatnonzero(self.listOfCheckboxes)
            self.selectedEpochs = self.epoch_var[self.nonZero]
            print(self.displayOption.get())
            self.displayOptions(self.displayOption)
            # self.plot_images(img)
        else:
            # there is no data loaded
            # dummy image; will never be displayed in canvas
            img = np.random.random((200,200))
            self.plot_images(img)



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
