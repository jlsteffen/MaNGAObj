"""
Created on Tue Sep 19 13:53:40 2017

This script will create a GUI which will allow the user to sort objects and 
non-objects into separate folders. This script is not strictly necessary, but
it makes the image sorting easier on the user.

True objects will have their png moved into the true directory and false 
objects will have their png moved into the false directory.

@author: Josh
"""

import tkinter as tk
import shutil as sy
import os, sys
from glob import glob

# Function to change image viewed in GUI
def next_img(event=' '):
    img_file = next(imgs)
    print(img_file)
    global img_name
    img_name = img_file
    img_label.img = tk.PhotoImage(file=img_file)
    img_label.config(image=img_label.img)
  
def f(event=' '):
    global img_name
    sy.move(str(img_name), str(false_dir))
    img_file = next(imgs)
    print(img_file)
    img_name = img_file
    img_label.img = tk.PhotoImage(file=img_file)
    img_label.config(image=img_label.img)

def t(event=' '):
    global img_name
    sy.move(str(img_name), str(true_dir))
    img_file = next(imgs)
    print(img_file)
    img_name = img_file
    img_label.img = tk.PhotoImage(file=img_file)
    img_label.config(image=img_label.img)

def close_window(event=' '): 
    root.destroy()
    
# Define directory with images and creates folders that the images will be 
# sorted into
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

true_dir = filepath + '/3.image_sort/true'
if not os.path.exists(true_dir):
    os.makedirs(true_dir)

false_dir = filepath + '/3.image_sort/false'
if not os.path.exists(false_dir):
    os.makedirs(false_dir)

os.chdir(filepath)
files = sorted(glob(filepath + '/3.image_sort/*.png'))
imgs = iter(files)

# Sets GUI
root = tk.Tk()

# Sets image and buttons in GUI
img_label = tk.Label(root)
img_label.pack(side='top')
frame = tk.Frame(root)
frame.pack()

# Move image to /false
button1 = tk.Button(frame, text='TRUE <A>', width=25, command=t)
button1.focus()
button1.pack()
# Move image to /true 
button2 = tk.Button(frame, text='FALSE <D>', width=25, command=f)
button2.pack()
next_img()
button3 = tk.Button(frame, text='QUIT <ESC>', width = 25, command=close_window)
button3.pack()

root.bind('<a>', t)
root.bind('<d>', f)
root.bind('<Escape>', close_window)

root.mainloop()