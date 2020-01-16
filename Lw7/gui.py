import l3 as l
import math as mth
import tkinter as tk
from tkinter import ttk
from tkinter import *
from PIL import Image, ImageTk
LARGE_FONT = ("Verdana", 12)


class Application(tk.Tk):

    def __init__(self, *args, **kwargs):  # конструктор класса
        # позиционные и именные аргументы
        tk.Tk.__init__(self, *args, **kwargs)

        label = tk.Label(self, text="Лабораторная №3", font=LARGE_FONT)
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

        label1 = tk.Label(self, text='Погрешность')
        label1.pack()
        self.ent1 = tk.Entry(self, width=20, bd=3)
        self.ent1.pack()

        label4 = tk.Label(self, text='Выберете метод')
        label4.pack()
        self.cb1 = ttk.Combobox(self,
                                values=[u"Метод простых итераций"
                                        " (метод Либмана)",
                                        u"Метод Зейделя",
                                        u"Метод простых итераций"
                                        " с верхней релаксацией"],
                                height=3)
        self.cb1.set('Не выбрано')
        self.cb1.pack()

        self.button = ttk.Button(
            self, text="Построить!",
            command=self.GetParams)

        self.button.pack(pady=10, padx=10)
        self.button.pack()
        #self.ent0.insert(0,3)
        # self.ent1.insert(0,1)
        self.ent2.insert(0,20)
        self.ent3.insert(0,20)
        self.ent1.insert(0,0.00001)

    def ErrorWindow(self, note):
        err_win = Tk()
        err_note = tk.Label(
            err_win, text=note,
            width=50, height=10)
        err_win.title('ERROR')
        err_note.pack(pady=10, padx=10)

    def GetParams(self):
        x0 = 0
        xl = mth.pi

        y0 = 0
        yl = 1

        x = int(self.ent2.get())
        y = int(self.ent3.get())
        method = self.cb1.get()
        eps = float(self.ent1.get())

        hx = (xl - x0) / (x)
        hy = (yl - y0) / (y)

        if method == 'Метод простых итераций (метод Либмана)':
            l.Make_Graph('simple', hx, hy, eps)
        elif method == "Метод Зейделя":
            l.Make_Graph('zeidel', hx, hy, eps)
        elif method == "Метод простых итераций с верхней релаксацией":
            l.Make_Graph('relax', hx, hy, eps)


if __name__ == '__main__':
    app = Application()
    app.title('Лабораторная №3. (c) Шевчук П.В.')
    app.mainloop()
