"""
Задача 1(а) — Линеаризация.

Метод:
  Пространство состояний  ΔV = [y, ẏ, θ, θ̇]  (отклонения от номинала)

       | 0  1   0     0 |        | 0 |
  A =  | 0  0  g_lat  0 |,  B = | 0 |
       | 0  0   0     1 |        | 0 |
       | 0  0   0     0 |        | 1 |

  Управление с задержкой:
      ΔU(t) = −K · ΔV(t−τ)
  где K = [0, 0, a_ctrl, 0] (PD по углу θ)

  Уравнение движения:
      ΔV' = A·ΔV + B·ΔU  →  θ'' + a·θ(t−τ) = 0
                           ÿ = g_lat·θ

  Интегрирование: явный метод Эйлера (1-й порядок).
  Задержка реализована через буфер истории состояний.

Начальные условия (φ₀ = 0 — наихудший случай по углу):
    θ(0) = A,  θ'(0) = 0,  y(0) = 0,  ẏ(0) = 0

Выход:
  Консоль: боковое смещение на высоте 10 км.
  CSV:     out_task1_linearization.csv  (t, theta, omega, y, vy, H)
"""

import math
import csv
import sys
import params as P


def run(A=P.A, tau=P.tau, T_rock=P.T_rock, h=P.h_step, a_ctrl=P.a_ctrl,
        phi0=0.0, out_csv='out_task1_linearization.csv'):
    """
    Запуск линеаризованной симуляции.

    Parameters
    ----------
    phi0 : float
        Начальная фаза качки (рад).
    out_csv : str
        Имя выходного CSV-файла.

    Returns
    -------
    float
        Боковое смещение y (м) на высоте H_10km.
    """
    Omega = 2 * math.pi / T_rock
    t_end = P.time_to_height(P.H_10km)
    n_steps = int(t_end / h) + 2
    n_delay = max(1, round(tau / h))   # число шагов задержки

    # Начальные условия  ΔV = [y, vy, theta, omega]
    state = [0.0,
             0.0,
             A * math.cos(phi0),
             -A * Omega * math.sin(phi0)]

    # Буфер истории θ (для доступа к θ(t−τ))
    theta_buf = [state[2]] * (n_delay + 1)

    # Матрица A (по строкам)
    g_lat = P.g_lat
    A_mat = [
        [0.0, 1.0, 0.0,    0.0],   # ẏ  = vy
        [0.0, 0.0, g_lat,  0.0],   # ÿ  = g_lat·θ
        [0.0, 0.0, 0.0,    1.0],   # θ' = ω
        [0.0, 0.0, 0.0,    0.0],   # ω' = ΔU
    ]
    B = [0.0, 0.0, 0.0, 1.0]      # ΔU воздействует на ω'
    K = [0.0, 0.0, a_ctrl, 0.0]   # ΔU = −K·ΔV(t−τ)

    rows = []  # для CSV (заполняется только если out_csv не None)
    y_at_10km = None
    t = 0.0

    for i in range(n_steps):
        H_now = P.height_at(t)
        if H_now >= P.H_10km and y_at_10km is None:
            y_at_10km = state[0]

        if out_csv is not None:
            rows.append((t, state[2], state[3], state[0], state[1], H_now))

        # Задержанное состояние (из буфера)
        delayed = theta_buf[0]

        # Управление: ΔU(t) = −K·ΔV(t−τ); используем только K[2]·θ(t−τ)
        delta_U = -(K[2] * delayed)

        # Производная  dΔV/dt = A·ΔV + B·ΔU  (явный Эйлер)
        deriv = [
            sum(A_mat[row][col] * state[col] for col in range(4)) + B[row] * delta_U
            for row in range(4)
        ]

        state = [state[j] + h * deriv[j] for j in range(4)]
        t += h

        # Обновить буфер истории (FIFO)
        theta_buf.pop(0)
        theta_buf.append(state[2])

    if y_at_10km is None:
        y_at_10km = state[0]

    # Запись CSV
    if out_csv is not None:
        with open(out_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['t_s', 'theta_rad', 'omega_rad_s', 'y_m', 'vy_m_s', 'H_m'])
            for row in rows:
                writer.writerow([f'{v:.6f}' for v in row])

    return y_at_10km


def main():
    args = P.parse_args('Задача 1(а): линеаризация, отклонение на 10 км')
    print('=' * 60)
    print("Задача 1(а): Линеаризация  ΔV' = A·ΔV + B·ΔU,  ΔU = −K·ΔV(t−τ)")
    print('=' * 60)
    print(f'  A      = {args.A} рад  ({math.degrees(args.A):.2f}°)')
    print(f'  τ      = {args.tau} с')
    print(f'  T_rock = {args.T} с')
    print(f'  h_step = {args.h} с')
    print(f'  a_ctrl = {args.a} рад/с²/рад')
    print(f'  Время до H=10 км: {P.time_to_height(P.H_10km):.2f} с')
    print()
    print('Матрица A пространства состояний  (ΔV = [y, ẏ, θ, θ̇]):')
    print('  | 0   1    0       0  |   (ẏ   = vy)')
    print('  | 0   0  g_lat     0  |   (ÿ   = g_lat·θ)')
    print("  | 0   0    0       1  |   (θ'  = ω)")
    print("  | 0   0    0       0  |   (ω'  = ΔU)")
    print(f'  g_lat = {P.g_lat} м/с²,  B = [0,0,0,1]ᵀ,  K = [0,0,{args.a},0]')
    print()

    # Sweep по фазам для полной картины
    print(f'  {"φ₀ (°)":>8}  {"y на 10 км (м)":>16}')
    print(f'  {"-"*8}  {"-"*16}')
    Omega = 2 * math.pi / args.T
    max_y = 0.0
    for k in range(8):
        phi0 = 2 * math.pi * k / 8
        csv_file = 'out_task1_linearization.csv' if k == 0 else None
        y = run(A=args.A, tau=args.tau, T_rock=args.T, h=args.h, a_ctrl=args.a,
                phi0=phi0, out_csv=csv_file)
        print(f'  {math.degrees(phi0):>8.1f}  {y:>16.2f}')
        max_y = max(max_y, abs(y))

    print(f'  {"Максимум |y|":>8}  {max_y:>16.2f}')
    print()
    print('  Временной ряд для φ₀=0° сохранён в out_task1_linearization.csv')


if __name__ == '__main__':
    main()
