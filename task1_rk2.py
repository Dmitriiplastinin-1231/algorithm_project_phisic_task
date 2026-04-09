"""
Задача 1(б) — RK2 с задержкой (метод Хойна / предиктор–корректор).

Модель:
    θ'' + a·θ(t−τ) = 0        (угловая динамика с задержанным управлением)
    ÿ   = g_lat · θ(t)        (боковое ускорение)

Начальные условия (φ₀ задаётся аргументом, по умолчанию φ₀=0):
    θ(0) = A·cos(φ₀),   θ'(0) = −A·Ω·sin(φ₀),   y(0)=0,  ẏ(0)=0

Метод Хойна (RK2):
    k₁ = f(tₙ, Vₙ)          — θ_задержан = θ(tₙ−τ) из буфера
    k₂ = f(tₙ+h, Vₙ+h·k₁) — θ_задержан = θ(tₙ+h−τ) = θ(tₙ−τ+h) из буфера
    Vₙ₊₁ = Vₙ + h/2·(k₁+k₂)

Буфер истории θ: кольцевой массив длиной (n_delay+2).
Задержанное значение при шаге i: theta_buf[i − n_delay].

Выход:
    Консоль: боковое смещение на H=10 км, сравнение с линеаризацией.
    CSV:     out_task1_rk2.csv  (t, theta, omega, y, vy, H)
"""

import math
import csv
import params as P


def _make_f(g_lat, a_ctrl):
    """Возвращает функцию правой части f(state, theta_delayed)."""
    def f(state, theta_del):
        y, vy, theta, omega = state
        return (vy, g_lat * theta, omega, -a_ctrl * theta_del)
    return f


def run(A=P.A, tau=P.tau, T_rock=P.T_rock, h=P.h_step, a_ctrl=P.a_ctrl,
        phi0=0.0, H_target=P.H_10km, out_csv='out_task1_rk2.csv'):
    """
    Запуск RK2-симуляции с задержкой.

    Parameters
    ----------
    phi0 : float
        Начальная фаза качки (рад).
    H_target : float
        Целевая высота (м); интеграция выполняется до достижения H_target.
    out_csv : str or None
        Имя CSV-файла; если None — файл не создаётся.

    Returns
    -------
    float
        Боковое смещение y (м) на высоте H_target.
    """
    Omega = 2 * math.pi / T_rock
    t_end = P.time_to_height(H_target)
    n_steps = int(t_end / h) + 2
    n_delay = max(1, round(tau / h))   # количество шагов задержки

    # Начальные условия  V = [y, vy, theta, omega]
    state = [0.0,
             0.0,
             A * math.cos(phi0),
             -A * Omega * math.sin(phi0)]

    # Буфер истории θ (длина n_delay+2 для индекса i-n_delay и i-n_delay+1)
    buf_size = n_delay + 2
    theta_buf = [state[2]] * buf_size   # до старта θ = θ(0)

    f = _make_f(P.g_lat, a_ctrl)

    rows = []
    y_at_target = None
    t = 0.0

    for i in range(n_steps):
        H_now = P.height_at(t)
        if H_now >= H_target and y_at_target is None:
            y_at_target = state[0]

        if out_csv is not None:
            rows.append((t, state[2], state[3], state[0], state[1], H_now))

        # Задержанные θ для шагов k₁ (tₙ) и k₂ (tₙ+h)
        buf_idx = i % buf_size
        # θ(tₙ−τ): позиция в кольцевом буфере на n_delay шагов назад
        idx_del_k1 = (buf_idx - n_delay) % buf_size
        idx_del_k2 = (buf_idx - n_delay + 1) % buf_size
        theta_del_k1 = theta_buf[idx_del_k1]
        theta_del_k2 = theta_buf[idx_del_k2]

        # ── Метод Хойна (RK2) ──────────────────────────────────────────
        k1 = f(state, theta_del_k1)
        state_pred = [state[j] + h * k1[j] for j in range(4)]
        k2 = f(state_pred, theta_del_k2)
        state = [state[j] + h * 0.5 * (k1[j] + k2[j]) for j in range(4)]
        t += h

        # Записать новый θ в кольцевой буфер
        theta_buf[(buf_idx + 1) % buf_size] = state[2]

    if y_at_target is None:
        y_at_target = state[0]

    # Запись CSV
    if out_csv is not None:
        with open(out_csv, 'w', newline='', encoding='utf-8') as fout:
            writer = csv.writer(fout)
            writer.writerow(['t_s', 'theta_rad', 'omega_rad_s', 'y_m', 'vy_m_s', 'H_m'])
            for row in rows:
                writer.writerow([f'{v:.6f}' for v in row])

    return y_at_target


def main():
    args = P.parse_args('Задача 1(б): RK2 с задержкой, отклонение на 10 км')
    print('=' * 60)
    print('Задача 1(б): RK2 (Хойн) с буфером задержки τ')
    print('=' * 60)
    print(f'  A      = {args.A} рад  ({math.degrees(args.A):.2f}°)')
    print(f'  τ      = {args.tau} с')
    print(f'  T_rock = {args.T} с')
    print(f'  h_step = {args.h} с')
    print(f'  a_ctrl = {args.a} рад/с²/рад')
    print(f'  φ₀     = 0 рад  (наихудший случай по θ)')
    print(f'  Время до H=10 км: {P.time_to_height(P.H_10km):.2f} с')
    print()

    y_rk2 = run(A=args.A, tau=args.tau, T_rock=args.T, h=args.h, a_ctrl=args.a,
                phi0=0.0, H_target=P.H_10km, out_csv='out_task1_rk2.csv')

    # Сравнение с линеаризацией (Эйлер)
    import task1_linearization as lin
    y_lin = lin.run(A=args.A, tau=args.tau, T_rock=args.T, h=args.h, a_ctrl=args.a,
                    phi0=0.0, out_csv=None)

    print(f'  RK2  (Хойн):         y = {y_rk2:.2f} м')
    print(f'  Линеаризация (Эйлер): y = {y_lin:.2f} м')
    if abs(y_lin) > 1e-9:
        diff_pct = abs(y_rk2 - y_lin) / abs(y_lin) * 100
        print(f'  Расхождение:          {diff_pct:.1f} %')
    print()
    print('  Алгоритм RK2 (Хойн):')
    print('    k₁ = f(tₙ, Vₙ,   θ(tₙ−τ))')
    print('    k₂ = f(tₙ+h, Vₙ+h·k₁, θ(tₙ+h−τ))')
    print('    Vₙ₊₁ = Vₙ + h/2·(k₁+k₂)')
    print()

    # Сравнение методов по всем фазам
    import task1_linearization as lin
    print(f'  {"φ₀ (°)":>8}  {"RK2 y (м)":>12}  {"Эйлер y (м)":>12}  {"разн. %":>9}')
    print(f'  {"-"*8}  {"-"*12}  {"-"*12}  {"-"*9}')
    max_rk2 = 0.0
    for k in range(8):
        phi0 = 2 * math.pi * k / 8
        yr = run(A=args.A, tau=args.tau, T_rock=args.T, h=args.h, a_ctrl=args.a,
                 phi0=phi0, H_target=P.H_10km, out_csv=None)
        ye = lin.run(A=args.A, tau=args.tau, T_rock=args.T, h=args.h, a_ctrl=args.a,
                     phi0=phi0, out_csv=None)
        dp = abs(yr - ye) / (abs(yr) + 1e-12) * 100
        print(f'  {math.degrees(phi0):>8.1f}  {yr:>12.2f}  {ye:>12.2f}  {dp:>9.1f}')
        max_rk2 = max(max_rk2, abs(yr))
    print(f'  {"Макс |y|":>8}  {max_rk2:>12.2f}')
    print()
    print('  Временной ряд для φ₀=0° сохранён в out_task1_rk2.csv')


if __name__ == '__main__':
    main()
