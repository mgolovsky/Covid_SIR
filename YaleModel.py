import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T_end = 501
Uninfected = np.zeros(T_end)
Uninfected[0] = 4990

Asympt = np.zeros(T_end)
Asympt[0] = 10

Fp = np.zeros(T_end)  # false pos
Tp = np.zeros(T_end)  # true pos
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
X = 5/21



def plot(U, A, S, R,D, Fp, Tp, t):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    # ax.plot(t, G, 'b-', label='Oral OP Concentration (μg/L)', markevery=10)
    #ax.plot(t, U, 'g-', label='Uninfected')
    ax.plot(t, A, 'b-', label='Asymptomatic')
    ax.plot(t, S, 'r-', label='Symptomatic')
    ax.plot(t, R, 'm-', label='Recovered')
    ax.plot(t, Fp, 'c-', label='False Pos')
    ax.plot(t, Tp, 'y-', label='True Pos')
    ax.plot(t, D, 'k-', label='Dead')


    # Set labels and turn grid on
    ax.set(xlabel='Time $t$, 8 hour intervals', ylabel=r'Population', title='Students in each Pop')
    ax.grid(True)
    ax.legend(loc='upper left')
    # Use space most effectively
    fig.tight_layout()
    # Save as a PNG file
    #fig.savefig('Oseltamivir_Concentration_{}.png'.format(plot_name))
    fig.show()

def InfectionModel(t):

    Uninfected[t + 1] = Uninfected[t] * (1 - Beta * (Asympt[t] / (Uninfected[t] + Asympt[t]))) \
                        - Uninfected[t - 1] * screening * (1 - specificity) + mu * Fp[t] - X
    Asympt[t + 1] = Asympt[t] * (
            (1 - symptom_onset - recovery_rate) + Beta * (Uninfected[t] / (Uninfected[t] + Asympt[t]))) \
                    - Asympt[t - 1] * screening * sensitivity + X

    Fp[t + 1] = Fp[t] * (1 - mu) + Uninfected[t - 1] * screening * (1 - specificity)

    Tp[t + 1] = Tp[t] * (1 - symptom_onset - recovery_rate) + Asympt[t - 1] * screening * sensitivity

    Sympt[t + 1] = Sympt[t] * (1 - recovery_rate - death) + symptom_onset * (Tp[t] + Asympt[t])

    Recovered[t + 1] = Recovered[t] + recovery_rate * (Tp[t] + Asympt[t] + Sympt[t])
    print(Tp[t])

    Deaths[t + 1] = Deaths[t] + death * Sympt[t]


time = range(0, T_end)
for t in time[0:500]:
    InfectionModel(t)

plot(Uninfected,Asympt,Sympt,Recovered,Deaths, Fp, Tp, time)

# out = odeint(InfectionModel, Yo, t)
