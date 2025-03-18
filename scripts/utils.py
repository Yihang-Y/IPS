import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML, display
from matplotlib.gridspec import GridSpec
import os
import tempfile
from PIL import Image
from tqdm import tqdm

def kinetic_energy(velocity):
    return 0.5 * np.linalg.norm(velocity) ** 2

def pot_func_LJ(dist, epsilon, sigma):
    return 4 * epsilon * ((sigma / dist) ** 12 - (sigma / dist) ** 6)

def potential_energy(position, epsilon=1.0, sigma=1.0, pot = "LJ"):
    num_particles = len(position[0])
    total_potential_energy = 0
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            pos_i = np.array([position[0][i], position[1][i]])
            pos_j = np.array([position[0][j], position[1][j]])
            dist = np.linalg.norm(pos_i - pos_j)
            if pot == "LJ":
                total_potential_energy += pot_func_LJ(dist, epsilon, sigma)
            else:
                raise ValueError("Potential not implemented")
    return total_potential_energy

def calculate_total_kinetic_energy(velocities):
    num_particles = len(velocities[0])
    total_kinetic_energy = 0
    for i in range(num_particles):
        vel = np.array([velocities[0][i], velocities[1][i]])
        total_kinetic_energy += kinetic_energy(vel)
    return total_kinetic_energy

def calculate_total_energy(positions, velocities):
    num_particles = len(positions[0])
    total_kinetic_energy = 0
    total_potential_energy = 0
    for i in range(num_particles):
        vel = np.array([velocities[0][i], velocities[1][i]])
        total_kinetic_energy += kinetic_energy(vel)
    total_potential_energy = potential_energy(positions)
    return total_kinetic_energy + total_potential_energy

class SimulationVisualizer:
    def __init__(self, simulator, particle_system, 
                 rad=10.0, draw_interval=100, dt=0.001):
        """
        visualizer for the simulation
        :param simulator: IPSimulator object
        :param particle_system: ParticleSystem object
        :param rad: constraint radius
        :param draw_interval: draw interval
        :param dt: time step
        """
        self.simulator = simulator
        self.p = particle_system
        self.rad = rad
        self.dt = dt
        self.draw_interval = draw_interval
        self.callbacks = []
        self.fig = plt.figure(figsize=(12, 6))
        
        
        self._init_energy_containers()
        
    def _init_plots(self):
        """init the plots"""
        
        self.ax1 = self.fig.add_subplot(121)
        self.scat = self.ax1.scatter([], [], s=50, edgecolors='k')
        self._draw_constraint()
        self.ax1.set_xlim(-self.rad-2, self.rad+2)
        self.ax1.set_ylim(-self.rad-2, self.rad+2)
        self.ax1.set_aspect('equal')
        self.ax1.set_title(f'Particle Positions (R={self.rad})')
        
        self.ax2 = self.fig.add_subplot(122)
        self.lines = {
            'total': self.ax2.plot([], [], 'k-', label='Total')[0],
            'kinetic': self.ax2.plot([], [], 'r--', label='Kinetic')[0],
            'potential': self.ax2.plot([], [], 'b-.', label='Potential')[0]
        }
        self.ax2.set_title('Energy Evolution')
        self.ax2.legend()
        return self.scat, *self.lines.values()

    def _draw_constraint(self):
        """draw the constraint circle"""
        constraint = plt.Circle((0,0), self.rad, color='r', fill=False, 
                               linestyle='--', alpha=0.7)
        self.ax1.add_artist(constraint)
        
    def _init_energy_containers(self):
        """generate containers for data"""
        self.energy_history = {
            'total': [],
            'kinetic': [],
            'potential': []
        }
        
    def add_callback(self, callback_func):
        """
        add user defined callback function
        :param callback_func: 
        """
        self.callbacks.append(callback_func)
        
    def _update_energy_plot_limits(self):
        """adjust the limits of the energy plot"""
        all_energies = np.concatenate([
            self.energy_history['total'],
            self.energy_history['kinetic'],
            self.energy_history['potential']
        ])
        
        if len(all_energies) > 0:
            pad = 0.1 * (np.max(all_energies) - np.min(all_energies))
            self.energy_ylim = [
                np.min(all_energies) - pad,
                np.max(all_energies) + pad
            ]
            self.ax2.set_ylim(*self.energy_ylim)
            
    def _calculate_energies(self):
        """calculate the energies"""
        positions = self.p.get_positions()
        velocities = self.p.get_velocities()
        
        total_energy = calculate_total_energy(positions, velocities)
        kinetic = total_energy - potential_energy(positions)
        potential = potential_energy(positions)
                
        return {
            'total': total_energy,
            'kinetic': kinetic,
            'potential': potential
        }
    
    def _update_plots(self, step):
        """update the plots"""
        # update particle positions
        positions = self.p.get_positions()
        self.scat.set_offsets(np.c_[positions[0], positions[1]])
        
        # update energy plot
        energies = self._calculate_energies()
        for key in self.energy_history:
            self.energy_history[key].append(energies[key])
            self.lines[key].set_data(
                range(len(self.energy_history[key])),
                self.energy_history[key]
            )
            
        # adjust limits
        self._update_energy_plot_limits()
        self.ax2.set_xlim(0, len(self.energy_history['total']))
        
    def run_animation(self, total_steps):
        # print(f"Running animation for {total_steps} steps")

        progress_bar = tqdm(total=total_steps // self.draw_interval, desc="Animating", unit="frame")

        """do the animation"""
        def _frame_update(frame):
            # do the callbacks
            for cb in self.callbacks:
                cb(self.simulator, frame*self.draw_interval)
                
            self._update_plots(frame)
            # do the simulation
            self.simulator.integrate_n_steps(self.dt, self.draw_interval)
            progress_bar.update(1)
            return self.scat, *self.lines.values()
        
        ani = animation.FuncAnimation(
            self.fig, _frame_update,
            frames=range(0, total_steps, self.draw_interval), init_func=self._init_plots,
        )
        json = ani.to_jshtml()
        progress_bar.close()
        plt.close()
        return HTML(json)
    
    def run_animation_to_gif(self, total_steps, gif_path, fps=10):
        """
        output the animation to a gif
        :param total_steps: total steps
        :param gif_path: output gif path
        :param fps: frame per second
        """
        print(f"Running animation for {total_steps} steps, saving to {gif_path}")
        
        # create a temporary directory to store the frames
        temp_dir = tempfile.mkdtemp()
        frame_files = []
        self._init_plots()
        try:
            # generate frames
            for frame_idx, step in enumerate(range(0, total_steps, self.draw_interval)):
                for cb in self.callbacks:
                    cb(self.simulator, step)

                self._update_plots(frame_idx)
                
                self.simulator.integrate_n_steps(self.dt, self.draw_interval)

                frame_path = os.path.join(temp_dir, f"frame_{frame_idx:04d}.png")
                self.fig.savefig(frame_path, dpi=100)
                frame_files.append(frame_path)

        finally:
            print(f"Output number of frames: {len(frame_files)}")
            # out put a gif
            self._create_gif(frame_files, gif_path, fps)
            for frame_file in frame_files:
                os.remove(frame_file)
            os.rmdir(temp_dir)

    def _create_gif(self, frame_files, gif_path, fps):
        
        images = []
        for frame_file in frame_files:
            img = Image.open(frame_file)
            images.append(img)

        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            optimize=True,
            duration=1000 // fps,  
            loop=0
        )
        print(f"GIF saved to {gif_path}")
