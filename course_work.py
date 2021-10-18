# Автор: Попов Даниил
# Группа 131-Мко
from math import factorial, sin, log2, pi
from cmath import cosh, sinh
from tkinter import *
from tkinter import ttk


# Функция произведения матриц (для поиска степеней матриц)
def multi_matrix(A, B):
    multi_AB = []
    for i in range(len(A)):
        multi_row = []
        for j in range(len(A)):
            multi_row.append(0)
        multi_AB.append(multi_row)
    for i in range(len(A)):
        for j in range(len(A)):
            for k in range(len(A)):
                multi_AB[i][j] += A[i][k] * B[k][j]
    return multi_AB


# Функция умножения матрицы на число
def multi_matrix_number(A, k):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] = A[i][j] * k
    return A


# Функция сложения матриц
def sum_matrix(A, B):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] += B[i][j]
    return A


# Функция нахождения коэффициентов p характеристического многочлена
def coefficients_p(A):
    p = []
    s = []
    first_degree_A = A
    for i in range(len(A)):
        s.append(trace_matrix(first_degree_A))
        first_degree_A = multi_matrix(first_degree_A, A)
        if i == 0:
            p.append(s[i])
        else:
            p.append(s[i])
            k = i - 1
            for j in range(i):
                p[i] -= p[j] * s[k]
                k -= 1
            p[i] /= (i + 1)
    return p


# Функция поиска следа матрицы
def trace_matrix(A):
    tr = 0
    for i in range(len(A)):
        tr += A[i][i]
    return tr


def B_in_kappa (g, n, p):
    if g >= 0 and g <= (n - 2):
        return 0
    elif g == (n - 1):
        return 1
    elif g >= n:
        b_sum = 0
        for l in range(1, n):
            b_sum += p[l - 1] * B_in_kappa(g - l, n, p)
        return b_sum
    else:
        return "ERROR"


def search_for_kappa (index, n, Np, p):
    kappa = 1 / factorial(index)
    for g in range(index):
        temp_p = p[n - index + g - 1]
        for j in range(n, (n + Np)):
            temp_sum = B_in_kappa(j - 1 - g, n, p) / factorial(j)
        kappa += temp_p * temp_sum
    return kappa


# Основная часть программы
# Численное решение
Np = 11
n = 4
# lamda_0 - длина волны
lambda_0 = 830
k0 = 2 * pi / lambda_0
theta = pi / 3
# d - толщина слоя
d = 500
G = complex(0.4322, 0.0058)
eps1 = eps2 = eps3 = complex(-4.8984, 19.415)
im = complex(0, 1)
Wd = [
    [0, (1 - (sin(theta) ** 2 / eps3)) * im * k0, 0, 0],
    [eps1 * im * k0, 0, k0 * G, 0],
    [0, 0, 0, im * k0],
    [-k0 * G, 0, im * k0 * (eps2 - sin(theta) ** 2), 0]
]

# Теорема
# Поиск максимального элемента матрицы Wd
maxWd = 0
for i in range(n):
    for j in range(n):
        if abs(Wd[i][j]) >= abs(maxWd):
            maxWd = Wd[i][j]

# Поиск m
m = 2

while True:
    temp_beta = abs(maxWd) * d * (n + 4) / 2 / m
    if temp_beta < 1:
        break
    else:
        m *= 2


def no_param(*args):
    try:
        m = int(param_m.get())
        while True:
            m *= 2
            if abs(maxWd) * (n + 4) / 2 / m < 1:
                param_m.set(m)
                break
    except ValueError:
        pass


# Полиномиальная аппроксимация
def yes_param(*args):
    A = []
    m = int(param_m.get())
    for i in range(n):
        temp_row = []
        for j in range(n):
            temp_row.append(Wd[i][j] / m)
        A.append(temp_row)

    p = coefficients_p(A)

    first_degree_A = A
    for i in range(n):
        if i == 0:
            identity_matrix = []
            for l in range(len(A)):
                identity_row = []
                for k in range(len(A)):
                    if l == k:
                        identity_row.append(1)
                    else:
                        identity_row.append(0)
                identity_matrix.append(identity_row)
            exp_A = multi_matrix_number(identity_matrix, search_for_kappa(i, n, Np, p))
        else:
            exp_A = sum_matrix(exp_A, \
                multi_matrix_number(A, search_for_kappa(i, n, Np, p)))
            A = multi_matrix(A, first_degree_A)
    # Возведение exp_A в степень m
    for i in range(int(log2(m))):
        exp_A = multi_matrix(exp_A, exp_A)
    exp_Wd = exp_A
    # Таблица 1
    for i in range(n):
        for j in range(n):
            e = Entry(mainframe, font=('Arial', 10))
            e.grid(row=i+3, column=j, sticky=(W, E))
            if exp_Wd[i][j].imag < 0:
                val = f"{'%.4f' % exp_Wd[i][j].real} - " \
                      f"{'%.4f' % abs(exp_Wd[i][j].imag)}i"
            else:
                val = f"{'%.4f' % exp_Wd[i][j].real} + " \
                      f"{'%.4f' % exp_Wd[i][j].imag}i"
            e.insert(END, val)


def exact(*args):
    # точное решение
    # z - количетсво слоев
    z = 1

    sigma_2 = k0 ** 2 * ((1 - ((sin(theta) ** 2) / eps3)) * eps1 +
                         (eps2 - (sin(theta) ** 2) - (G * G) / eps3))
    sigma_4 = k0 ** 4 * eps1 * ((G * G + eps2 * sin(theta) ** 2 \
            - sin(theta) ** 4) / eps3 - eps2 + sin(theta) ** 2)

    beta_1 = (sigma_2 / 2 + (sigma_2 ** 2 / 4 + sigma_4) ** 0.5) ** 0.5
    beta_2 = (sigma_2 / 2 - (sigma_2 ** 2 / 4 + sigma_4) ** 0.5) ** 0.5

    Sigma_0 = (beta_1 ** 2 * cosh(z * beta_2) - beta_2 ** 2
               * cosh(z * beta_1)) / (beta_1 ** 2 - beta_2 ** 2)
    Sigma_1 = (beta_1 ** 3 * sinh(z * beta_2) - beta_2 ** 3
               * sinh(z * beta_1)) / (beta_1 * beta_2 *
                                      (beta_1 ** 2 - beta_2 ** 2))
    Sigma_2 = (cosh(z * beta_2) - cosh(z * beta_1)) \
              / (beta_1 ** 2 - beta_2 ** 2)
    Sigma_3 = (beta_1 * sinh(z * beta_2) - beta_2 *
               sinh(z * beta_1)) / (beta_1 * beta_2 *
                                    (beta_1 ** 2 - beta_2 ** 2))

    w2 = multi_matrix(Wd, Wd)
    w3 = multi_matrix(w2, Wd)

    t11 = Sigma_0 + Sigma_1 * Wd[0][0] + Sigma_2 \
          * w2[0][0] + Sigma_3 * w3[0][0]
    t12 = Sigma_1 * Wd[0][1] + Sigma_2 * w2[0][1] + Sigma_3 * w3[0][1]
    t13 = Sigma_1 * Wd[0][2] + Sigma_2 * w2[0][2] + Sigma_3 * w3[0][2]
    t14 = Sigma_1 * Wd[0][3] + Sigma_2 * w2[0][3] + Sigma_3 * w3[0][3]
    t21 = Sigma_1 * Wd[1][0] + Sigma_2 * w2[1][0] + Sigma_3 * w3[1][0]
    t22 = Sigma_0 + Sigma_1 * Wd[1][1] + Sigma_2 \
          * w2[1][1] + Sigma_3 * w3[1][1]
    t23 = Sigma_1 * Wd[1][2] + Sigma_2 * w2[1][2] + Sigma_3 * w3[1][2]
    t24 = Sigma_1 * Wd[1][3] + Sigma_2 * w2[1][3] + Sigma_3 * w3[1][3]
    t31 = Sigma_1 * Wd[2][0] + Sigma_2 * w2[2][0] + Sigma_3 * w3[2][0]
    t32 = Sigma_1 * Wd[2][1] + Sigma_2 * w2[2][1] + Sigma_3 * w3[2][1]
    t33 = Sigma_0 + Sigma_1 * Wd[2][2] + Sigma_2 \
          * w2[2][2] + Sigma_3 * w3[2][2]
    t34 = Sigma_1 * Wd[2][3] + Sigma_2 * w2[2][3] + Sigma_3 * w3[2][3]
    t41 = Sigma_1 * Wd[3][0] + Sigma_2 * w2[3][0] + Sigma_3 * w3[3][0]
    t42 = Sigma_1 * Wd[3][1] + Sigma_2 * w2[3][1] + Sigma_3 * w3[3][1]
    t43 = Sigma_1 * Wd[3][2] + Sigma_2 * w2[3][2] + Sigma_3 * w3[3][2]
    t44 = Sigma_0 + Sigma_1 * Wd[3][3] + Sigma_2\
          * w2[3][3] + Sigma_3 * w3[3][3]

    T = [[t11, t12, t13, t14],
         [t21, t22, t23, t24],
         [t31, t32, t33, t34],
         [t41, t42, t43, t44]]
    # Таблица 2
    for i in range(n):
        for j in range(n):
            e = Entry(mainframe, font=('Arial', 10))
            e.grid(row=i + 9, column=j, sticky=(W, E))
            if T[i][j].imag < 0:
                val = f"{'%.4f' % T[i][j].real} - " \
                      f"{'%.4f' % abs(T[i][j].imag)}i"
            else:
                val = f"{'%.4f' % T[i][j].real} + " \
                      f"{'%.4f' % T[i][j].imag}i"
            e.insert(END, val)


root = Tk()
root.geometry("590x360")
root.title('Course Work. Popov Daniil.')
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

numeric_title = "Матрица переноса (численное вычисление)"
param_m_text = "Устроит ли параметр масштабирования m ="
param_m = StringVar()
param_m.set(m)
exact_title = "Матрица переноса (точное решение)"
style = ttk.Style()
style.map("TButton",
          foreground=[('pressed', 'black'), ('active', 'blue')],
          background=[('pressed', '!disabled', 'black'),
                      ('active', 'white')]
          )

ttk.Label(mainframe, text=numeric_title, font=('arial', 18),
          foreground='black', background="white",
          relief='solid').grid(row=0, column=0, columnspan=4, sticky=W)
ttk.Label(mainframe, text=param_m_text, font=('arial', 14),
          foreground='black').grid(column=0, row=1, columnspan=3, sticky=W)
ttk.Label(mainframe, textvariable=param_m, font=('arial', 14),
          foreground='black', relief='groove', width=5)\
    .grid(column=3, row=1, sticky=W)
ttk.Button(mainframe, text="Да", command=yes_param)\
    .grid(column=0, row=2, sticky=W)
ttk.Button(mainframe, text="Нет", command=no_param)\
    .grid(column=1, row=2, sticky=W)
ttk.Label(mainframe, text=exact_title, font=('arial', 18),
          foreground='black', background="white",
          relief='solid').grid(row=7, column=0, columnspan=4,
                               pady=50, sticky=W)
ttk.Button(mainframe, text="Вычислить", command=exact)\
    .grid(column=0, row=8, sticky=W)

for child in mainframe.winfo_children():
    child.grid_configure(padx=5, pady=5)

root.mainloop()
