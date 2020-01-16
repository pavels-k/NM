import l4 as l
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import *
from PIL import Image, ImageTk
LARGE_FONT = ("Verdana", 12)


class Application(tk.Tk):

    def __init__(self, *args, **kwargs):  # конструктор класса
        # позиционные и именные аргументы
        tk.Tk.__init__(self, *args, **kwargs)

        label = tk.Label(self, text="Лабораторная №4", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        im = Image.open('img.gif')
        render = ImageTk.PhotoImage(im)
        img = tk.Label(self, image=render)
        img.image = render
        img.pack(pady=10, padx=10)

        label2 = tk.Label(self, text='Введите число интервалов по X')
        label2.pack()
        self.ent2 = tk.Entry(self, width=20, bd=3)
        self.ent2.pack()

        label3 = tk.Label(self, text='Введите число интервалов по Y')
        label3.pack()
        self.ent3 = tk.Entry(self, width=20, bd=3)
        self.ent3.pack()

        label31 = tk.Label(self, text='Введите число интервалов по T')
        label31.pack()
        self.ent31 = tk.Entry(self, width=20, bd=3)
        self.ent31.pack()

        label0 = tk.Label(self, text='Параметр a > 0')
        label0.pack()
        self.ent0 = tk.Entry(self, width=20, bd=3)
        self.ent0.pack()

        label1 = tk.Label(self, text='Параметр t')
        label1.pack()
        self.ent1 = tk.Entry(self, width=20, bd=3)
        self.ent1.pack()

        label4 = tk.Label(self, text='Выберете метод')
        label4.pack()
        self.cb1 = ttk.Combobox(self,
                                values=[u"Метод переменных направлений",
                                        u"Метод дробных шагов"],
                                height=3)
        self.cb1.set('Не выбрано')
        self.cb1.pack()

        self.button = ttk.Button(
            self, text="Построить!",
            command=self.GetParams)

        self.button.pack(pady=10, padx=10)
        self.button.pack()

        self.ent0.insert(0,3)
        self.ent2.insert(0,20)
        self.ent3.insert(0,20)
        self.ent31.insert(0,20)
        self.ent1.insert(0,0.5)

    def ErrorWindow(self, note):
        err_win = Tk()
        err_note = tk.Label(
            err_win, text=note,
            width=50, height=10)
        err_win.title('ERROR')
        err_note.pack(pady=10, padx=10)

    def GetParams(self):
        x0 = 0
        xl = np.pi / 4

        y0 = 0
        yl = np.log(2)

        x_int = int(self.ent2.get())
        y_int = int(self.ent3.get())
        r_int = int(self.ent31.get())
        par_a = float(self.ent0.get())
        t = float(self.ent1.get())
        method = self.cb1.get()

        if x_int < y_int:
            tmp = y_int
            y_int = x_int
            x_int = tmp
        hx = (xl - x0) / (x_int)
        hy = (yl - y0) / (y_int)
        r = (t) / (r_int)

        if method == 'Метод переменных направлений':
            l.Make_Graph("Alternating", hx, hy, r, t, par_a)
        elif method == "Метод дробных шагов":
            l.Make_Graph("Fractional", hx, hy, r, t, par_a)


if __name__ == '__main__':
    app = Application()
    app.title('Лабораторная №4. (c) Шевчук П.В.')
    app.mainloop()
