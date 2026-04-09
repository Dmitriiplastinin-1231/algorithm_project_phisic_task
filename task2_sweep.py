"""
Задача 2 — Максимальное боковое смещение при выведении на орбиту (160 км).

Метод:
  RK2 (Хойн) с буфером задержки + перебор фаз φ₀ ∈ [0, 2π).

Начальные условия для каждой фазы:
    θ(0) = A·cos(φ₀),   θ'(0) = −A·Ω·sin(φ₀)

Sweep:
  N_small = 8   — инженерная оценка (равномерный шаг 45°)
  N_large = 72  — верификационный перебор (шаг 5°)
  Сравниваются полученные максимумы.

Выход:
  Консоль: таблица [φ₀, y_max_N8, y_max_N72], итоговый max.
  CSV:     out_task2_sweep_N{N}.csv  (phi0_deg, y_at_orbit_m)
           out_task2_comparison.csv  (N, max_deviation_m)
"""

import math
import csv
import params as P
from task1_rk2 import run as rk2_run


def sweep(N, A=P.A, tau=P.tau, T_rock=P.T_rock, h=P.h_step, a_ctrl=P.a_ctrl,
          H_target=P.H_160km):
    """
    Перебор N равномерно распределённых фаз φ₀ ∈ [0, 2π).

    Returns
    -------
    list of (phi0_deg, y_at_target)
    """
    results = []
    for k in range(N):
        phi0 = 2 * math.pi * k / N
        y = rk2_run(A=A, tau=tau, T_rock=T_rock, h=h, a_ctrl=a_ctrl,
                    phi0=phi0, H_target=H_target, out_csv=None)
        results.append((math.degrees(phi0), y))
    return results


def max_abs_deviation(results):
    """Максимальное |y| из списка результатов."""
    return max(abs(y) for _, y in results)


def save_csv(results, filename):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['phi0_deg', 'y_at_orbit_m'])
        for phi0_deg, y in results:
            writer.writerow([f'{phi0_deg:.4f}', f'{y:.4f}'])


def main():
    args = P.parse_args('Задача 2: sweep фаз + RK2, max смещение на 160 км')
    print('=' * 60)
    print('Задача 2: Sweep фаз φ₀, RK2 с задержкой, H = 160 км')
    print('=' * 60)
    print(f'  A      = {args.A} рад  ({math.degrees(args.A):.2f}°)')
    print(f'  τ      = {args.tau} с')
    print(f'  T_rock = {args.T} с')
    print(f'  h_step = {args.h} с')
    print(f'  a_ctrl = {args.a} рад/с²/рад')
    print(f'  Время до H=160 км: {P.time_to_height(P.H_160km):.2f} с')
    print()

    N_small = P.N_phases_small
    N_large = P.N_phases_large

    # ── Базовый перебор N=8 ───────────────────────────────────────────────
    print(f'  Запуск sweep N={N_small} фаз (шаг {360//N_small}°)...')
    res_small = sweep(N_small, A=args.A, tau=args.tau, T_rock=args.T,
                      h=args.h, a_ctrl=args.a)
    max_small = max_abs_deviation(res_small)
    save_csv(res_small, f'out_task2_sweep_N{N_small}.csv')
    print(f'  {'φ₀ (°)':>10}  {'y на 160 км (м)':>18}')
    print(f'  {"-"*10}  {"-"*18}')
    for phi_d, y in res_small:
        print(f'  {phi_d:>10.1f}  {y:>18.2f}')
    print(f'  Максимум |y| (N={N_small}): {max_small:.2f} м')
    print()

    # ── Расширенный перебор N=72 ──────────────────────────────────────────
    print(f'  Запуск sweep N={N_large} фаз (шаг {360//N_large}°)...')
    res_large = sweep(N_large, A=args.A, tau=args.tau, T_rock=args.T,
                      h=args.h, a_ctrl=args.a)
    max_large = max_abs_deviation(res_large)
    save_csv(res_large, f'out_task2_sweep_N{N_large}.csv')
    print(f'  Максимум |y| (N={N_large}): {max_large:.2f} м')
    print()

    # ── Сравнение N=8 и N=72 ─────────────────────────────────────────────
    print('  Сравнение N=8 и N=72:')
    print(f'  {"N":>5}  {"max |y| (м)":>14}  {"Примечание"}')
    print(f'  {"-"*5}  {"-"*14}  {"-"*30}')
    print(f'  {N_small:>5}  {max_small:>14.2f}  инженерная оценка')
    print(f'  {N_large:>5}  {max_large:>14.2f}  верификация')
    if max_large > 1e-9:
        err_pct = abs(max_small - max_large) / max_large * 100
        print(f'  Расхождение: {err_pct:.1f}%')
    print()

    # ── Итоговый ответ ────────────────────────────────────────────────────
    max_dev = max(max_small, max_large)
    print(f'  *** Расчётный максимум бокового смещения: {max_dev:.2f} м ***')
    if max_dev < 1000.0:
        print('  → Смещение < 1 км: ракета укладывается в допуск ±1 км')
    else:
        print('  → Смещение > 1 км: превышение допуска ±1 км!')

    # ── CSV для сравнения ─────────────────────────────────────────────────
    with open('out_task2_comparison.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['N_phases', 'max_deviation_m'])
        writer.writerow([N_small, f'{max_small:.4f}'])
        writer.writerow([N_large, f'{max_large:.4f}'])

    print()
    print('  Результаты сохранены:')
    print(f'    out_task2_sweep_N{N_small}.csv')
    print(f'    out_task2_sweep_N{N_large}.csv')
    print('    out_task2_comparison.csv')


if __name__ == '__main__':
    main()
