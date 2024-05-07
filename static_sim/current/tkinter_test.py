# Python program to create 
# a file explorer in Tkinter
  
# import all components
# from the tkinter library
from tkinter import *
  
# import filedialog module
from tkinter import filedialog, ttk
from tkinter.messagebox import showinfo
import os
                                                                                             
# Create the root window
window = Tk()
# Set window title
window.title('Electrode Simulation')
# Set window size
# window.geometry()

'''
Trap type selecton.
'''
def show_selected_trap_type():
    showinfo(
        title='Result',
        message=selected_trap_type.get()
    )
    print(selected_trap_type.get())
    if(selected_trap_type.get() == trap_types[1][1]):
        label_file_explorer.pack(fill='x', padx=5, pady=5)
        button_explore.pack(fill='x', padx=5, pady=5)

trap_types = (('Approximate Quetzal (Symmetric).', 'Approx'),
         ('GDS', 'GDS'))

# label
label = ttk.Label(text="Select trap type:")
label.pack(fill='x', padx=5, pady=5)

# set default selected trap type.
selected_trap_type = StringVar(value=trap_types[0][1])
# radio buttons
def radio_command():
    if(selected_trap_type.get() == trap_types[1][1]):
        label_file_explorer.pack(fill='x', padx=5, pady=5)
        button_explore.pack(fill='x', padx=5, pady=5)
    else:
        label_file_explorer.pack_forget()
        button_explore.pack_forget()

for trap_type in trap_types:
    r = ttk.Radiobutton(
        window,
        text=trap_type[0],
        value=trap_type[1],
        variable=selected_trap_type,
        command=radio_command
    )
    r.pack(fill='x', padx=5, pady=5)

# button
# button = ttk.Button(
#     window,
#     text="Get Selected Trap",
#     command=show_selected_trap_type)
# button.pack(fill='x', padx=5, pady=5)
  
#Set window background color
# window.config(background = "gray")
  
# Function for opening the 
# file explorer window
gds_filename = ''
def browseFiles():
    global gds_filename
    gds_filename = filedialog.askopenfilename(initialdir = os.getcwd(),
                                          title = 'Select a File',
                                          filetypes = (('GDS files',
                                                        '*.gds'),
                                                       ('All files',
                                                        '*.*')))
      
    # Change label contents
    label_file_explorer.configure(text="File Opened: " + gds_filename)
# Create a File Explorer label

# if (selected_trap_type.get() == trap_types[1][1]):
label_file_explorer = Label(window, 
                            text = "Select trap GDS:",
                            height = 4, 
                            fg = "white")
# label_file_explorer.pack(fill='x', padx=5, pady=5)
    
button_explore = ttk.Button(window, 
                        text = "Browse Files",
                        command = browseFiles) 
# button_explore.pack(fill='x', padx=5, pady=5)

# button_print = Button(window, 
#                     text = "Print filename",
#                     command = lambda: print(gds_filename)) 

# button_exit = Button(window, 
#                      text = "Exit",
#                      command = exit) 
# button_exit.pack(fill='x', padx=5, pady=5)

# Grid method is chosen for placing
# the widgets at respective positions 
# in a table like structure by
# specifying rows and columns
# label_file_explorer.grid(column = 1, row = 1)
  
# button_explore.grid(column = 2, row = 1)
  
# button_exit.grid(column = 1,row = 3)
  
# Let the window wait for any events
window.mainloop()