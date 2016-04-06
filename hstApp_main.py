#!/usr/bin/env python

import Tkinter as tk
import os


class TopLevelWin:

    def __init__(self, master):
        self.master = master

        self.label1 = tk.Label(self.master, text='What would you like to do?')
        self.label1.pack()

        self.frame = tk.Frame(self.master)

        self.button1 = tk.Button(self.frame, text='Analyze radial profiles', width=25,
                                 command=lambda: os.system('python hstApp_mod1.py'))
        self.button1.pack()

        self.button2 = tk.Button(self.frame, text='Analyze iso-velocity images', width=25,
                                 command=lambda: os.system('python hstApp_mod2.py'))
        self.button2.pack()
        #
        # self.button3 = tk.Button(self.frame, text='Analyze E0', width=25,
        #                          command=lambda: os.system('python AP_E0.py'))
        # self.button3.pack()
        #
        # self.button4 = tk.Button(self.frame, text='specER', width=25)
        #
        # self.button4.pack()

        self.quitbutton = tk.Button(self.frame, text='Quit hstApp', width=25,
                                    command=self.frame.quit)
        self.quitbutton.pack()

        self.frame.pack()


def main():
    root = tk.Tk()
    root.title('hstApp - Main module')
    root.resizable(width=tk.FALSE, height=tk.FALSE)
    # root.minsize(width=200, height=200)
    # root.maxsize(width=200, height=200)
    TopLevelWin(root)
    root.mainloop()


if __name__ == '__main__':
    __author__ = 'M. Teodoro'
    main()
