# Jain, Kritika
# 1001-093-381
# 2017-09-24
# Assignment_02_02


import tkinter as tk
from tkinter import simpledialog
from tkinter import filedialog


class cl_widgets:
    def __init__(self, ob_root_window, ob_world=[]):
        self.ob_root_window = ob_root_window
        self.ob_world = ob_world
        self.buttons_panel_01 = cl_buttons_panel_01(self)
        self.buttons_panel_02 = cl_buttons_panel_02(self)
        self.statusBar_frame = cl_statusBar_frame(self.ob_root_window)
        self.statusBar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.statusBar_frame.set('%s', 'This is the status bar')
        self.ob_canvas_frame = cl_canvas_frame(self)
        self.ob_world.add_canvas(self.ob_canvas_frame.canvas)


class cl_canvas_frame:
    def __init__(self, master):
        self.master = master
        self.canvas = tk.Canvas(master.ob_root_window, width=880, height=880, bg="yellow")
        self.canvas.pack(expand=tk.YES, fill=tk.BOTH)
        self.canvas.bind('<Configure>', self.canvas_resized_callback)
        self.canvas.bind("<ButtonPress-1>", self.left_mouse_click_callback)
        self.canvas.bind("<ButtonRelease-1>", self.left_mouse_release_callback)
        self.canvas.bind("<B1-Motion>", self.left_mouse_down_motion_callback)
        self.canvas.bind("<ButtonPress-3>", self.right_mouse_click_callback)
        self.canvas.bind("<ButtonRelease-3>", self.right_mouse_release_callback)
        self.canvas.bind("<B3-Motion>", self.right_mouse_down_motion_callback)
        self.canvas.bind("<Key>", self.key_pressed_callback)
        self.canvas.bind("<Up>", self.up_arrow_pressed_callback)
        self.canvas.bind("<Down>", self.down_arrow_pressed_callback)
        self.canvas.bind("<Right>", self.right_arrow_pressed_callback)
        self.canvas.bind("<Left>", self.left_arrow_pressed_callback)
        self.canvas.bind("<Shift-Up>", self.shift_up_arrow_pressed_callback)
        self.canvas.bind("<Shift-Down>", self.shift_down_arrow_pressed_callback)
        self.canvas.bind("<Shift-Right>", self.shift_right_arrow_pressed_callback)
        self.canvas.bind("<Shift-Left>", self.shift_left_arrow_pressed_callback)
        self.canvas.bind("f", self.f_key_pressed_callback)
        self.canvas.bind("b", self.b_key_pressed_callback)

    def key_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Key pressed')

    def up_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Up arrow was pressed")

    def down_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Down arrow was pressed")

    def right_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Right arrow was pressed")

    def left_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Left arrow was pressed")

    def shift_up_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Shift up arrow was pressed")

    def shift_down_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Shift down arrow was pressed")

    def shift_right_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Shift right arrow was pressed")

    def shift_left_arrow_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "Shift left arrow was pressed")

    def f_key_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "f key was pressed")

    def b_key_pressed_callback(self, event):
        self.master.statusBar_frame.set('%s', "b key was pressed")

    def left_mouse_click_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Left mouse button was clicked. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = event.x
        self.y = event.y
        self.canvas.focus_set()

    def left_mouse_release_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Left mouse button was released. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = None
        self.y = None

    def left_mouse_down_motion_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Left mouse down motion. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = event.x
        self.y = event.y

    def right_mouse_click_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Right mouse down motion. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = event.x
        self.y = event.y

    def right_mouse_release_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Right mouse button was released. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = None
        self.y = None

    def right_mouse_down_motion_callback(self, event):
        self.master.statusBar_frame.set('%s', 'Right mouse down motion. ' + \
                                        'x=' + str(event.x) + '   y=' + str(event.y))
        self.x = event.x
        self.y = event.y

    def canvas_resized_callback(self, event):

        self.master.statusBar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.master.statusBar_frame.set('%s', 'Canvas width = ' + str(self.canvas.cget("width")) + \
                                        '   Canvas height = ' + str(self.canvas.cget("height")))
        self.canvas.pack()
        self.master.ob_world.redisplay(self.master.ob_canvas_frame.canvas)


class cl_buttons_panel_01:
    def __init__(self, master):
        self.master = master
        frame = tk.Frame(master.ob_root_window)
        frame.pack()
        self.label_file = tk.Label(frame, text="File")
        self.label_file.pack(side=tk.LEFT)

        self.file = tk.Entry(frame)
        self.file.pack(side=tk.LEFT)


        self.file_dialog_button = tk.Button(frame, text="Open", fg="black", command=self.load)
        self.file_dialog_button.pack(side=tk.LEFT)

        self.line = tk.Label(frame, text="Line")
        self.line.pack(side=tk.LEFT)


        self.start_point = tk.Entry(frame)
        self.start_point.insert(0, "0,0,0")
        self.start_point.pack(side=tk.LEFT)

        self.label_end_point = tk.Label(frame, text=",")
        self.label_end_point.pack(side=tk.LEFT)

        self.end_point = tk.Entry(frame)
        self.end_point.insert(0, "0,0,1")
        self.end_point.pack(side=tk.LEFT)

        self.label_deg = tk.Label(frame, text="deg")
        self.label_deg.pack(side=tk.LEFT)


        self.input_degrees = tk.Entry(frame)
        self.input_degrees.insert(0, "90")
        self.input_degrees.pack(side=tk.LEFT)

        self.label_steps = tk.Label(frame, text="steps")
        self.label_steps.pack(side=tk.LEFT)

        self.line_input_Steps = tk.Entry(frame)
        self.line_input_Steps.insert(0, "100")
        self.line_input_Steps.pack(side=tk.LEFT)

        self.rotate = tk.Button(frame, text="rotate", command=self.rotate_func)
        self.rotate.pack(side=tk.LEFT)

        self.var_filename = tk.StringVar()
        self.var_filename.set('')




    def rotate_func(self):  # this function prepares data for rotation
        angle = int(self.input_degrees.get()) / int(self.line_input_Steps.get())
        star_point = self.start_point.get()
        star_point_list = star_point.split(",")
        end_point = self.end_point.get()
        end_point_list = end_point.split(",")
        for i in range(1,int(self.line_input_Steps.get())+1):
            self.master.ob_world.vertexes = self.master.ob_world.rotate(angle, star_point_list, end_point_list, self.master.ob_world.vertexes).copy()
            self.master.ob_canvas_frame.canvas_resized_callback(None)



    def say_hi(self):
        self.master.statusBar_frame.set('%s', "hi there, everyone!")

    def ask_for_string(self):
        s = simpledialog.askstring('My Dialog', 'Please enter a string')
        self.master.statusBar_frame.set('%s', s)

    def ask_for_float(self):
        f = float(simpledialog.askfloat('My Dialog', 'Please enter a float'))

        self.master.statusBar_frame.set('%s', str(f))

    def browse_file(self):
        self.var_filename.set(tk.filedialog.askopenfilename(filetypes=[("allfiles", "*"), ("pythonfiles", "*.txt")]))
        filename = self.var_filename.get()
        self.master.statusBar_frame.set('%s', filename)


    def load(self):  # data is loaded from file in here
        self.var_filename.set(tk.filedialog.askopenfilename(filetypes=[("allfiles", "*"), ("pythonfiles", "*.txt")]))
        filename = self.var_filename.get()
        self.master.statusBar_frame.set('%s', filename)
        self.file.insert(0, filename)
        vertexes = []
        faces = []
        window = []
        viewport = []
        with open(filename) as f:
            file_lines = f.readlines()
        file_lines = [x.strip() for x in file_lines]
        for line in file_lines:
            if line != "" :  # line by line analysis
                if line[0] == "f":
                    faces.append([int(x) for x in cl_buttons_panel_01.prepare_data(line)])
                if line[0] == "w":
                    window = cl_buttons_panel_01.prepare_data(line)
                if line[0] == "s":
                    viewport = cl_buttons_panel_01.prepare_data(line)
                if line[0] == 'v':
                    vertexes.append(cl_buttons_panel_01.prepare_data(line))
        self.master.ob_world.vertexes = vertexes
        self.master.ob_world.faces = faces        # data is translated up to the master object
        self.master.ob_world.window = window
        self.master.ob_world.viewport = viewport
        self.draw() # initial draw is made

    def draw(self):
        self.master.ob_canvas_frame.canvas.delete("all") # canvas cleared
        self.master.ob_world.create_graphic_objects(self.master.ob_canvas_frame.canvas)


    def prepare_data(str):
        a_list = str.split()
        a_list.pop(0)
        for i in a_list:
            a_list[a_list.index(i)] = float(i)
        return a_list




class cl_buttons_panel_02:
    def __init__(self, master):
        self.master = master
        frame = tk.Frame(master.ob_root_window)
        frame.pack()

        self.label1 = tk.Label(frame, text="point")
        self.label1.pack(side=tk.LEFT)

        self.input_point = tk.Entry(frame)  #text field for the input
        self.input_point.insert(0, "0,0,0")
        self.input_point.pack(side=tk.LEFT)

        self.label2 = tk.Label(frame, text="scale")
        self.label2.pack(side=tk.LEFT)

        self.input_scale = tk.Entry(frame)
        self.input_scale.insert(0, "1.1,1.1,1.1")
        self.input_scale.pack(side=tk.LEFT)


        self.label3 = tk.Label(frame, text="steps")
        self.label3.pack(side=tk.LEFT)

        self.input_steps1 = tk.Entry(frame)
        self.input_steps1.insert(0, "100")
        self.input_steps1.pack(side=tk.LEFT)

        self.rotate = tk.Button(frame, text="Scale", command=self.scale_callback)
        self.rotate.pack(side=tk.LEFT)
        self.label4 = tk.Label(frame, text="translate")
        self.label4.pack(side=tk.LEFT)

        self.input_translate = tk.Entry(frame)
        self.input_translate.insert(0, "1,1,1")
        self.input_translate.pack(side=tk.LEFT)

        self.label3 = tk.Label(frame, text="steps")
        self.label3.pack(side=tk.LEFT)

        self.input_steps = tk.Entry(frame)
        self.input_steps.insert(0, "100")
        self.input_steps.pack(side=tk.LEFT)

        self.rotate = tk.Button(frame, text="translate", command=self.translate)
        self.rotate.pack(side=tk.LEFT)

    def translate(self):  # translation is made here
        step = int(self.input_steps.get())  # input from step
        translate = self.input_translate.get()  # input from move
        translate = translate.split(",")
        translate = [float(x) / step for x in translate]   # translation is divided in to smaller steps

        for i in range(1, int(self.input_steps.get()) + 1):
            self.master.ob_world.vertexes = self.master.ob_world.translate_vector(translate,
                                                                                  self.master.ob_world.vertexes).copy()
            self.master.ob_canvas_frame.canvas_resized_callback(None)


    def scale_callback(self):  # scaling data colected here
        step = int(self.input_steps1.get())
        scale = self.input_scale.get()
        scale = scale.split(",")
        scale = [float(x)**(1/step) for x in scale]  # scale is divided into smaller steps as a root of original

        point = self.input_point.get()
        point = point.split(",")
        point = [float(x) / step for x in point]
        for i in range(1, int(self.input_steps1.get()) + 1):
            self.master.ob_world.vertexes = self.master.ob_world.scale(point,scale,self.master.ob_world.vertexes).copy()
            self.master.ob_canvas_frame.canvas_resized_callback(None)



    def open_dialog_callback(self):
        d = MyDialog(self.master.ob_root_window)
        self.master.statusBar_frame.set('%s', "mydialog_callback pressed. Returned results: " + str(d.result))





# class MyDialog(tk.simpledialog.Dialog):
#     def body(self, master):
#
#         tk.Label(master, text="Integer:").grid(row=0, sticky=tk.W)
#         tk.Label(master, text="Float:").grid(row=1, column=0, sticky=tk.W)
#         tk.Label(master, text="String:").grid(row=1, column=2, sticky=tk.W)
#         self.e1 = tk.Entry(master)
#         self.e1.insert(0, 0)
#         self.e2 = tk.Entry(master)
#         self.e2.insert(0, 4.2)
#         self.e3 = tk.Entry(master)
#         self.e3.insert(0, 'Default text')
#
#         self.e1.grid(row=0, column=1)
#         self.e2.grid(row=1, column=1)
#         self.e3.grid(row=1, column=3)
#
#         self.cb = tk.Checkbutton(master, text="Hardcopy")
#         self.cb.grid(row=3, columnspan=2, sticky=tk.W)
#
#     def apply(self):
#         try:
#             first = int(self.e1.get())
#             second = float(self.e2.get())
#             third = self.e3.get()
#             self.result = first, second, third
#         except ValueError:
#             tk.tkMessageBox.showwarning(
#                 "Bad input",
#                 "Illegal values, please try again"
#             )
#
#
class cl_statusBar_frame(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.label = tk.Label(self, bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.label.pack(fill=tk.X)

    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()


# class cl_menu:
#     def __init__(self, master):
#         self.master = master
#         self.menu = tk.Menu(master.ob_root_window)
#         master.ob_root_window.config(menu=self.menu)
#         self.filemenu = tk.Menu(self.menu)
#         self.menu.add_cascade(label="File", menu=self.filemenu)
#         self.filemenu.add_command(label="New", command=self.menu_callback)
#         self.filemenu.add_command(label="Open...", command=self.menu_callback)
#         self.filemenu.add_separator()
#         self.filemenu.add_command(label="Exit", command=self.menu_callback)
#         self.dummymenu = tk.Menu(self.menu)
#         self.menu.add_cascade(label="Dummy", menu=self.dummymenu)
#         self.dummymenu.add_command(label="Item1", command=self.menu_item1_callback)
#         self.dummymenu.add_command(label="Item2", command=self.menu_item2_callback)
#         self.helpmenu = tk.Menu(self.menu)
#         self.menu.add_cascade(label="Help", menu=self.helpmenu)
#         self.helpmenu.add_command(label="About...", command=self.menu_help_callback)
#
#     def menu_callback(self):
#         self.master.statusBar_frame.set('%s', "called the menu callback!")
#
#     def menu_help_callback(self):
#         self.master.statusBar_frame.set('%s', "called the help menu callback!")
#
#     def menu_item1_callback(self):
#         self.master.statusBar_frame.set('%s', "called item1 callback!")
#
#     def menu_item2_callback(self):
#         self.master.statusBar_frame.set('%s', "called item2 callback!")
#

