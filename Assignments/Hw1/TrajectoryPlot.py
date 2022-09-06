import pandas as pd
import matplotlib.pyplot as pyp

Traj = pd.read_csv("Trajectory.csv")
Traj.sort_values('time', inplace = True)

fig, (ax_h, ax_v) = pyp.subplots(2, sharex = True, dpi = 120)
ax_h.plot(Traj['time'], Traj['altitude'])
ax_h.set_ylabel('Altitude [m]')
ax_h.grid(True)

ax_v.plot(Traj['time'], Traj['velocity'])
ax_v.set_xlabel('Time [s]')
ax_v.set_ylabel('Velocity [m/s]')
ax_v.grid(True)

fig.suptitle("Rocket Trajectory")

fig.show()

print(Traj)