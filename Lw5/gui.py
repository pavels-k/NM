import l1 as l
import math as mth
import tkinter as tk
from tkinter import ttk
from tkinter import *
from PIL import Image, ImageTk
LARGE_FONT = ("Verdana", 12)


class Application(tk.Tk):

    def __init__(self, *args, **kwargs):  # конструктор класса

        tk.Tk.__init__(self, *args, **kwargs)

        label = tk.Label(self, text="Лабораторная №1", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        im = Image.open('img.gif')
        render = ImageTk.PhotoImage(im)
        img = tk.Label(self, image=render)
        img.image = render
        img.pack(pady=10, padx=10)

        label0 = tk.Label(self, text='Введите параметр "t" > 0')
        label0.pack()
        self.ent0 = tk.Entry(self, width=20, bd=3)
        self.ent0.pack()

        label1 = tk.Label(self, text='Введите параметр "a" > 0')
        label1.pack()
        self.ent1 = tk.Entry(self, width=20, bd=3)
        self.ent1.pack()

        label2 = tk.Label(self, text='Введите число интервалов по X')
        label2.pack()
        self.ent2 = tk.Entry(self, width=20, bd=3)
        self.ent2.pack()

        label3 = tk.Label(self, text='Введите число интервалов по T')
        label3.pack()
        self.ent3 = tk.Entry(self, width=20, bd=3)
        self.ent3.pack()

        label4 = tk.Label(self, text='Выберете метод')
        label4.pack()
        self.cb1 = ttk.Combobox(self,
                                values=[u"Явный",
                                        u"Неявный",
                                        u"Явный-Неявный"],
                                height=3)
        self.cb1.set('Не выбрано')
        self.cb1.pack()

        label5 = tk.Label(
            self, text='Выберете аппроксимацию граничных условий')
        label5.pack()
        self.cb2 = ttk.Combobox(self,
                                values=[u"Двухточечный(первый порядок)",
                                        u"Трехточечный(второй порядок)",
                                        u"Двухточечный(второй порядок)"],
                                height=3)
        self.cb2.set('Не выбрано')
        self.cb2.pack()

        label6 = tk.Label(self, text='Введите момент времени от 0 до t')
        label6.pack()
        self.ent4 = tk.Entry(self, width=20, bd=3)
        self.ent4.pack()
        self.ent0.insert(0,2)
        self.ent1.insert(0,1)
        self.ent2.insert(0,21)
        self.ent3.insert(0,2001)
        self.ent4.insert(0,0.01)
        self.button = ttk.Button(
            self, text="Построить!",
            command=self.GetParams)

        self.button.pack(pady=10, padx=10)
        self.button.pack()
        

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

        param_t = float(self.ent0.get())
        param_a = float(self.ent1.get())
        m = int(self.ent2.get())
        n = int(self.ent3.get())
        method = self.cb1.get()
        apr = self.cb2.get()
        ans_t = float(self.ent4.get())
        aprox = 0
        

        h = (xl - x0) / (m - 1)
        tau = param_t / (n - 1)
        if param_a < 0:
            note = 'Ошибка!\nПараметр "а" должен быть больше 0!'
            self.ErrorWindow(note)
            return
        elif param_t < 0:
            note = 'Ошибка!\nПараметр "t" должен быть больше 0!'
            self.ErrorWindow(note)
            return
        elif ans_t > param_t: 
            note = 'Ошибка!\nПараметр момент времени должен быть меньше "t"!'
            self.ErrorWindow(note)
            return
        elif param_a**2 * tau / h**2 > 0.5\
                and method == 'Явный':
            note = 'Ошибка!\nПри таких мараметрах Явный метод не устойчив!\n\
            Пожалуйста, измените параметры сетки.'
            self.ErrorWindow(note)
            return
        if apr == 'Двухточечный(первый порядок)':
            aprox = 1
        elif apr == 'Трехточечный(второй порядок)':
            aprox = 2
        elif apr == 'Двухточечный(второй порядок)':
            aprox = 3
        ans_t = 0
        param_a = 1
        if method == 'Явный':
            l.explicit(param_t, param_a, h,
                       tau, m, n, aprox, ans_t)
        elif method == "Неявный":
            l.implicit(param_t, param_a, h,
                       tau, m, n, aprox, ans_t)
        elif method == "Явный-Неявный":
            l.KN(param_t, param_a, h,
                 tau, m, n, aprox, ans_t)

app = Application()
app.title('Лабораторная №1. (c) Шевчук П.В.')
app.mainloop()
