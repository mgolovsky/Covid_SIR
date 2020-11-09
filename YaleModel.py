import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T_end = 500
Susceptible = np.zeros(T_end)
Susceptible[0] = 4990

Undetected = np.zeros(T_end)
Undetected[0] = 10

Iso_un = np.zeros(T_end)  # false pos
Iso_in = np.zeros(T_end)  # true pos
Sympt = np.zeros(T_end)  # Symptomatic
Recovered = np.zeros(T_end)
Deaths = np.zeros(T_end)
Shock_campus = [x * 5 for x in range(0, T_end) if x % 168 == 0]
Shock_off = [x * 5 for x in range(0, T_end) if x % 168 == 0]

init_uninfect = 4900

Beta = .085
screening = (2 / 21)
death = 0.00004
recovery_rate = (1 / 21)
symptom_onset = 0.0102
mu = (1 / 3)  # false positives returned
sensitivity = .8
specificity = .98
X = 2/21



def plot(U, A, S, R,D, Fp, Tp, t):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    # ax.plot(t, G, 'b-', label='Oral OP Concentration (μg/L)', markevery=10)
    #ax.plot(t, U, 'g-', label='Uninfected')
    t = np.divide(t,3)

    ax.plot(t, A, 'b-', label='Asymptomatic')
    ax.plot(t, S, 'r-', label='Symptomatic')
    ax.plot(t, R, 'm-', label='Recovered')
    ax.plot(t, Fp, 'c-', label='False Pos')
    ax.plot(t, Tp, 'y-', label='True Pos')
    ax.plot(t, D, 'k-', label='Dead')


    # Set labels and turn grid on
    ax.set(xlabel='Time $d$, Days', ylabel=r'Population', title='Students in each Pop')
    ax.grid(True)
    ax.legend(loc='upper left')
    # Use space most effectively
    fig.tight_layout()
    # Save as a PNG file
    #fig.savefig('Oseltamivir_Concentration_{}.png'.format(plot_name))
    fig.show()

def InfectionModel(t):

    Susceptible[t + 1] = Susceptible[t] * (1 - Beta * (Undetected[t] / (Susceptible[t] + Undetected[t]))) \
                         - Susceptible[t - 1] * screening * (1 - specificity) + mu * Iso_un[t] - X
    Undetected[t + 1] = Undetected[t] * (
            (1 - symptom_onset - recovery_rate) + Beta * (Susceptible[t] / (Susceptible[t] + Undetected[t]))) \
                        - Undetected[t - 1] * screening * sensitivity + X

    Iso_un[t + 1] = Iso_un[t] * (1 - mu) + Susceptible[t - 1] * screening * (1 - specificity)

    Iso_in[t + 1] = Iso_in[t] * (1 - symptom_onset - recovery_rate) + Undetected[t - 1] * screening * sensitivity

    Sympt[t + 1] = Sympt[t] * (1 - recovery_rate - death) + symptom_onset * (Iso_in[t] + Undetected[t])

    Recovered[t + 1] = Recovered[t] + recovery_rate * (Iso_in[t] + Undetected[t] + Sympt[t])

    Deaths[t + 1] = Deaths[t] + death * Sympt[t]

    if t % 21 == 0:
        print("True Pos:{}, Asympt: {}, Sympt: {}".format(Iso_in[t], Undetected[t], Sympt[t]))



time = range(0, T_end)

for t in time[0:T_end-1]:
    InfectionModel(t)

plot(Susceptible, Undetected, Sympt, Recovered, Deaths, Iso_un, Iso_in, time)

# out = odeint(InfectionModel, Yo, t)
