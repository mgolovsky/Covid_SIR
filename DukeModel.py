import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T_end = 100

# ON CAMPUS VARIABLES
Susceptible_c = np.zeros(T_end)
Susceptible_c[0] = 2495

Undetected_c = np.zeros(T_end)
Undetected_c[0] = 5

Iso_un_c = np.zeros(T_end)  # false pos
Iso_in_c = np.zeros(T_end)  # true pos
Sympt_c = np.zeros(T_end)  # Symptomatic
Recovered_c = np.zeros(T_end)
Deaths_c = np.zeros(T_end)

# OFF CAMPUS VARIABLES
Susceptible_o = np.zeros(T_end)
Susceptible_o[0] = 2495

Undetected_o = np.zeros(T_end)
Undetected_o[0] = 5

Iso_un_o = np.zeros(T_end)  # false pos
Iso_in_o = np.zeros(T_end)  # true pos
Sympt_o = np.zeros(T_end)  # Symptomatic
Recovered_o = np.zeros(T_end)
Deaths_o = np.zeros(T_end)

totSusceptible = np.zeros(T_end)
totUndetected = np.zeros(T_end)

# SWEEP VARIABLES
Beta_stu = .306  # Base Case   Each student to each student
Beta_c = .1  # Base Case  Between on campus students
screening = 0.2857  # Base Case
X = 2 / 7  # Base Case
Y = 4 / 7

init_uninfect = 4900

death = 0.00004
recovery_rate = (1 / 7)
symptom_onset = 0.06122
mu = (1 / 14)  # false positives returned
sensitivity = .8
specificity = .98


def plot(U, A, S, R, D, Fp, Tp, t, name):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    # ax.plot(t, G, 'b-', label='Oral OP Concentration (μg/L)', markevery=10)
    # ax.plot(t, U, 'g-', label='Uninfected')

    ax.plot(t, A, 'b-', label='Asymptomatic')
    ax.plot(t, S, 'r-', label='Symptomatic')
    ax.plot(t, R, 'm-', label='Recovered')
    ax.plot(t, Fp, 'c-', label='False Pos')
    ax.plot(t, Tp, 'y-', label='True Pos')
    ax.plot(t, D, 'k-', label='Dead')

    # Set labels and turn grid on
    ax.set(xlabel='Time $d$, Days', ylabel=r'Population', title=name)
    ax.grid(True)
    ax.legend(loc='upper left')
    # Use space most effectively
    fig.tight_layout()
    fig.savefig("{}.png".format(name.replace(' ', '_')))
    fig.show()


def OnCampus(t):
    infectedFrac = (
            totUndetected[t] / (totSusceptible[t] + totUndetected[t]))  # Fraction of non-isolated students infected

    infectedFrac_c = Undetected_c[t] / (Susceptible_c[t] + Undetected_c[t])

    print(infectedFrac, infectedFrac_c)

    Susceptible_c[t + 1] = Susceptible_c[t] * (1 - Beta_stu * infectedFrac - Beta_c * infectedFrac_c) \
                           - Susceptible_c[t - 1] * screening * (1 - specificity) + mu * Iso_un_c[t] - X

    Undetected_c[t + 1] = Undetected_c[t] * (1 - symptom_onset - recovery_rate) \
                          + Susceptible_c[t] *(Beta_stu * infectedFrac + Beta_c * infectedFrac_c)\
                          - Undetected_c[t - 1] * screening * sensitivity + X

    Iso_un_c[t + 1] = Iso_un_c[t] * (1 - mu) + Susceptible_c[t - 1] * screening * (1 - specificity)

    Iso_in_c[t + 1] = Iso_in_c[t] * (1 - symptom_onset - recovery_rate) + Undetected_c[t - 1] * screening * sensitivity

    Sympt_c[t + 1] = Sympt_c[t] * (1 - recovery_rate - death) + symptom_onset * (Iso_in_c[t] + Undetected_c[t])

    Recovered_c[t + 1] = Recovered_c[t] + recovery_rate * (Iso_in_c[t] + Undetected_c[t] + Sympt_c[t])

    Deaths_c[t + 1] = Deaths_c[t] + death * Sympt_c[t]

    # if t % 21 == 0:
    # print("True Pos:{}, Asympt: {}, Sympt: {}".format(Iso_in_c[t], Undetected_c[t], Sympt_c[t]))


def OffCampus(t):
    infectedFrac = (
            totUndetected[t] / (totSusceptible[t] + totUndetected[t]))  # Fraction of non-isolated students infected

    Susceptible_o[t + 1] = Susceptible_o[t] * (
            1 - Beta_stu * infectedFrac) - Susceptible_o[t - 1] * screening * (1 - specificity) + \
                           mu * Iso_un_o[t] - X

    Undetected_o[t + 1] = Undetected_o[t] * (1 - symptom_onset - recovery_rate) + \
                          Beta_stu * Susceptible_o[t] * infectedFrac \
                          - Undetected_o[t - 1] * screening * sensitivity + Y

    Iso_un_o[t + 1] = Iso_un_o[t] * (1 - mu) + Susceptible_o[t - 1] * screening * (1 - specificity)

    Iso_in_o[t + 1] = Iso_in_o[t] * (1 - symptom_onset - recovery_rate) + Undetected_o[t - 1] * screening * sensitivity

    Sympt_o[t + 1] = Sympt_o[t] * (1 - recovery_rate - death) + symptom_onset * (Iso_in_o[t] + Undetected_o[t])

    Recovered_o[t + 1] = Recovered_o[t] + recovery_rate * (Iso_in_o[t] + Undetected_o[t] + Sympt_o[t])

    Deaths_o[t + 1] = Deaths_o[t] + death * Sympt_o[t]

    # if t % 21 == 0:
    # print("True Pos:{}, Asympt: {}, Sympt: {}".format(Iso_in_o[t], Undetected_o[t], Sympt_o[t]))


time = range(0, T_end)

for t in time[0:T_end - 1]:
    totSusceptible[t] = Susceptible_c[t] + Susceptible_o[t]
    totUndetected[t] = Undetected_c[t] + Undetected_o[t]
    OnCampus(t)
    OffCampus(t)
    #print("Susceptible: {}\n Undetected: {}".format(totSusceptible[t], totUndetected[t]))

plot(Susceptible_c, Undetected_c, Sympt_c, Recovered_c, Deaths_c, Iso_un_c, Iso_in_c, time, "On Campus")
plot(Susceptible_o, Undetected_o, Sympt_o, Recovered_o, Deaths_o, Iso_un_o, Iso_in_o, time, "Off Campus")

# out = odeint(InfectionModel, Yo, t)
