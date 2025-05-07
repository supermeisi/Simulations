import numpy as np
import matplotlib.pyplot as plt
import math
import functions as f

# Parameters for fused silica and muon
X0_fused_silica = 12.55  # radiation length in cm
mass_muon = 105.66  # MeV/c^2
z_muon = 1  # muon charge
N = 100  # numbr of steps

# Step sizes from 10 µm to 1 mm (in cm)
dx_values_cm = np.linspace(0, 4, N)

# Particle momentum in GeV/c between 0.5 and 1.5 GeV/c
p_values_GeV = np.linspace(0.5, 1.5, N)
p_values_MeV = p_values_GeV * 1000

# Calculate theta values for muon in fused silica
theta_values_mrad = [[
    1 / math.sqrt(2) *
    f.highland_theta(p, f.beta(mass_muon, p), z_muon, dx, X0_fused_silica)
    for p in p_values_MeV
] for dx in dx_values_cm]

# Convert step sizes to mm for plotting/output
dx_values_mm = dx_values_cm * 10

# Plot
fig = plt.figure(figsize=(6, 4))
plt.imshow(theta_values_mrad,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[0.5, 1.5, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling.png', dpi=600)
plt.savefig('straggling.pdf')

theta_sigma_corrected = np.zeros((N, N))
theta_sigma_uncorrected = np.zeros((N, N))

L = 4  # length of the material in cm
p = 500  # momentum in GeV/c
beta = f.beta(mass_muon, p)

N_dx = 1000
dx = L / N_dx  # step size in cm

theta_history = []
theta_history_corrected = []

theta = 0.0
sigma = f.highland_theta(p, beta, z_muon, dx, X0_fused_silica)

np.random.seed(2)

for _ in range(N_dx):
    theta = np.random.normal(theta, sigma)  # random angular kick
    theta_history.append(theta)

theta_start = theta_history[0]
theta_end = theta_history[-1]
track_fit = np.linspace(theta_start, theta_end, N_dx)

theta_history = np.array(theta_history)

theta_history_corrected = theta_history - track_fit

# Plot
fig = plt.figure(figsize=(7, 4))
plt.plot(np.linspace(0, L, N_dx), theta_history, label='Uncorrected')
plt.plot(np.linspace(0, L, N_dx), theta_history_corrected, label='Corrected')
plt.xlabel("$x$ [cm]")
plt.ylabel("$\\theta$ [mrad]")
plt.title("Accumulated Angular Straggling")
plt.grid()
plt.tight_layout()
plt.savefig('straggling_step.png', dpi=600)
plt.savefig('straggling_step.pdf')
plt.close()

N_dx = 100  # number of steps
N_tracks = 100  # number of tracks to average over

with open('straggling.txt', 'w') as file:
    file.write("p L theta theta_track\n")
    # Loop over momentum values
    for i, p in enumerate(p_values_MeV):
        # Loop over step sizes
        for j, L in enumerate(dx_values_cm):
            id = i * N + j

            dx = L / N_dx  # step size in cm

            theta_tot = 0.0
            theta_tot_corrected = 0.0

            theta_history_seeds = []
            theta_history_seeds_corrected = []

            sigma = f.highland_theta(p, f.beta(mass_muon, p), z_muon, dx,
                                     X0_fused_silica)

            # Loop over different seeds
            for k in range(N_tracks):
                np.random.seed(k)

                # Initial direction angle (in radians)
                theta = 0.0
                theta_history = [theta]
                theta_history_corrected = [theta]

                # Step through the material
                for _ in range(N_dx):
                    theta = np.random.normal(theta,
                                             sigma)  # random angular kick
                    theta_history.append(theta)

                theta_start = theta_history[0]
                theta_end = theta_history[-1]
                track_fit = np.linspace(theta_start, theta_end, N_dx + 1)

                theta_history_corrected = theta_history - track_fit

                theta_history_seeds.append(theta_history)
                theta_history_seeds_corrected.append(theta_history_corrected)

                # Convert angle history to mrad
                theta_history_mrad = np.array(theta_history)
                theta_history_mrad_corrected = np.array(
                    theta_history_corrected)

                # print(theta_history_mrad.std(), theta_history_mrad_corrected.std())

                theta_tot += np.sqrt(np.mean(theta_history_mrad**2))
                theta_tot_corrected += np.sqrt(
                    np.mean(theta_history_mrad_corrected**2))

            theta_tot /= N_tracks
            theta_tot_corrected /= N_tracks

            theta_sigma_uncorrected[j, i] = theta_tot
            theta_sigma_corrected[j, i] = theta_tot_corrected

            print(p, L, theta_tot, theta_tot_corrected,
                  theta_tot_corrected / theta_tot)

            #file.write(f"{p} {L} {theta_tot} {theta_tot_corrected}\n")

            # print(theta_tot, theta_tot_corrected)

fig = plt.figure(figsize=(5, 4))
plt.imshow(theta_sigma_uncorrected,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[0.5, 1.5, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling_muon_uncorrected_step.png', dpi=600)
plt.savefig('straggling_muon_uncorrected_step.pdf')
plt.close()

fig = plt.figure(figsize=(6, 4))
plt.imshow(theta_sigma_corrected,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[0.5, 1.5, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling_muon_corrected_step.png', dpi=600)
plt.savefig('straggling_muon_corrected_step.pdf')
plt.close()

mass_pion = 139.57  # MeV/c^2
z_pion = 1  # pion charge

# Particle momentum in GeV/c between 0.5 and 1.5 GeV/c
p_values_GeV = np.linspace(1.5, 6.0, N)
p_values_MeV = p_values_GeV * 1000

with open('straggling.txt', 'w') as file:
    file.write("p L theta theta_track\n")
    # Loop over momentum values
    for i, p in enumerate(p_values_MeV):
        # Loop over step sizes
        for j, L in enumerate(dx_values_cm):
            id = i * N + j

            dx = L / N_dx  # step size in cm

            theta_tot = 0.0
            theta_tot_corrected = 0.0

            theta_history_seeds = []
            theta_history_seeds_corrected = []

            sigma = f.highland_theta(p, f.beta(mass_pion, p), z_pion, dx,
                                     X0_fused_silica)

            # Loop over different seeds
            for k in range(N_tracks):
                np.random.seed(k)

                # Initial direction angle (in radians)
                theta = 0.0
                theta_history = [theta]
                theta_history_corrected = [theta]

                # Step through the material
                for _ in range(N_dx):
                    theta = np.random.normal(theta,
                                             sigma)  # random angular kick
                    theta_history.append(theta)

                theta_start = theta_history[0]
                theta_end = theta_history[-1]
                track_fit = np.linspace(theta_start, theta_end, N_dx + 1)

                theta_history_corrected = theta_history - track_fit

                theta_history_seeds.append(theta_history)
                theta_history_seeds_corrected.append(theta_history_corrected)

                # Convert angle history to mrad
                theta_history_mrad = np.array(theta_history)
                theta_history_mrad_corrected = np.array(
                    theta_history_corrected)

                # print(theta_history_mrad.std(), theta_history_mrad_corrected.std())

                theta_tot += np.sqrt(np.mean(theta_history_mrad**2))
                theta_tot_corrected += np.sqrt(
                    np.mean(theta_history_mrad_corrected**2))

            theta_tot /= N_tracks
            theta_tot_corrected /= N_tracks

            theta_sigma_uncorrected[j, i] = theta_tot
            theta_sigma_corrected[j, i] = theta_tot_corrected

            print(p, L, theta_tot, theta_tot_corrected,
                  theta_tot_corrected / theta_tot)

            #file.write(f"{p} {L} {theta_tot} {theta_tot_corrected}\n")

            # print(theta_tot, theta_tot_corrected)

fig = plt.figure(figsize=(5, 4))
plt.imshow(theta_sigma_uncorrected,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[1.5, 6.0, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling_pion_uncorrected_step.png', dpi=600)
plt.savefig('straggling_pion_uncorrected_step.pdf')
plt.close()

fig = plt.figure(figsize=(6, 4))
plt.imshow(theta_sigma_corrected,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[1.5, 6.0, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling_pion_corrected_step.png', dpi=600)
plt.savefig('straggling_pion_corrected_step.pdf')
plt.close()
