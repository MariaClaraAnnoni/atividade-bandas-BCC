import numpy as np
import matplotlib.pyplot as plt

def place_labels(ax, x, values, texts, min_dy=0.15):
    pairs = sorted(zip(values, texts), key=lambda v: v[0])
    adjusted = []

    for val, text in pairs:
        y = val
        for prev_y, _ in adjusted:
            if abs(y - prev_y) < min_dy:
                y = prev_y + min_dy
        adjusted.append((y, text))
    
    for y, text in adjusted:
        ax.text(x, y, text, fontsize=7, ha='center', va='center',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0))

gamma_h = [
    [0,0,0], [0,-1,0], [1,0,0], [1,-1,-1], 
    [0,1,0], [0,-1,-1], [1,1,1], [-1,1,1]
]

h_p = [
    [0,0,0], [-1,0,0], [1,0,0], [1,-1,-1], 
    [0,1,0], [0,-1,-1], [-1,-1,0], [2,-1,-1], [-2,0,0]
]

p_gamma = [
    [0,0,0], [0,-1,0], [-2,0,1], [1,-1,-1], 
    [0,1,0], [0,-1,-1], [1,1,-1], [1,-1,0], [-2,0,0]
]

gamma_n = [
    [0,0,0], [0,-1,0], [1,0,0], [1,-1,-1], 
    [0,0,1], [0,-1,-1], [1,1,-1], [-1,1,1], 
    [1,-1,0], [0,0,-1], [-1,-1,0], [0,0,-2]
]

fig, ax = plt.subplots(figsize=(12, 8))
alfa = np.linspace(0, 1, 100)

cmap = plt.cm.get_cmap("tab20")

def plot_gamma_h(list, ax, alfa):
    mid_idx = len(alfa) // 2
    label_values = []
    label_texts = []

    for i, m in enumerate(list):
        m1, m2, m3 = m
        e = (alfa + m1)**2 + (m2)**2 + (m3)**2
        ax.plot(alfa, e, color=cmap(i % 20), linewidth=1.2)

        label_values.append(round(e[mid_idx], 2))
        label_texts.append(str(m))

    place_labels(ax, alfa[mid_idx], label_values, label_texts)


def plot_h_p(list, ax, alfa):
    mid_idx = len(alfa) // 2
    label_values = []
    label_texts = []

    for i, m in enumerate(list):
        m1, m2, m3 = m
        
        kx = 1 - 0.5 * alfa
        ky = 0.5 * alfa
        kz = 0.5 * alfa

        e = (kx + m1)**2 + (ky + m2)**2 + (kz + m3)**2
        alfa_plot = alfa + 1

        ax.plot(alfa_plot, e, color=cmap(i % 20), linewidth=1.2)

        label_values.append(round(e[mid_idx], 2))
        label_texts.append(str(m))

    place_labels(ax, alfa[mid_idx] + 1, label_values, label_texts)


def plot_p_gamma(list, ax, alfa):
    mid_idx = len(alfa) // 2
    label_values = []
    label_texts = []

    for i, m in enumerate(list):
        m1, m2, m3 = m
        
        k_comp = 0.5 - 0.5 * alfa
        e = (k_comp + m1)**2 + (k_comp + m2)**2 + (k_comp + m3)**2

        alfa_plot = alfa + 2
        ax.plot(alfa_plot, e, color=cmap(i % 20), linewidth=1.2)

        label_values.append(round(e[mid_idx], 2))
        label_texts.append(str(m))

    place_labels(ax, alfa[mid_idx] + 2, label_values, label_texts)


def plot_gamma_n(list, ax, alfa):
    mid_idx = len(alfa) // 2
    label_values = []
    label_texts = []

    for i, m in enumerate(list):
        m1, m2, m3 = m
        
        kx = 0.5 * alfa
        ky = 0.5 * alfa
        kz = 0

        e = (kx + m1)**2 + (ky + m2)**2 + (kz + m3)**2

        alfa_plot = alfa + 3
        ax.plot(alfa_plot, e, color=cmap(i % 20), linewidth=1.2)

        label_values.append(round(e[mid_idx], 2))
        label_texts.append(str(m))

    place_labels(ax, alfa[mid_idx] + 3, label_values, label_texts)


plot_gamma_h(gamma_h, ax, alfa)
plot_h_p(h_p, ax, alfa)
plot_p_gamma(p_gamma, ax, alfa)
plot_gamma_n(gamma_n, ax, alfa)
ax.set_ylim(0, 5)
ax.set_xlim(0, 4)

ax.axvline(x=1, color='black', linestyle='--', linewidth=0.5)
ax.axvline(x=2, color='black', linestyle='--', linewidth=0.5)
ax.axvline(x=3, color='black', linestyle='--', linewidth=0.5)

ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xticklabels([r'$\Gamma$', r'$H$', r'$P$', r'$\Gamma$', r'$N$'])

ax.set_ylabel(r'Energia Normalizada $\epsilon / [\frac{\hbar^2}{2m}(\frac{2\pi}{a})^2]$', fontsize=12)
ax.set_title('Estrutura de Bandas BCC - Aproximação de Rede Vazia')
ax.grid(True, linestyle=':', alpha=0.5)

plt.show()
