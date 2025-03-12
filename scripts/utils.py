import numpy as np

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

def calculate_total_energy(positions, velocities):
    num_particles = len(positions[0])
    total_kinetic_energy = 0
    total_potential_energy = 0
    for i in range(num_particles):
        vel = np.array([velocities[0][i], velocities[1][i]])
        total_kinetic_energy += kinetic_energy(vel)
    total_potential_energy = potential_energy(positions)
    return total_kinetic_energy + total_potential_energy

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML, display
from matplotlib.gridspec import GridSpec
import os
import tempfile
from PIL import Image

class SimulationVisualizer:
    def __init__(self, simulator, particle_system, 
                 rad=10.0, draw_interval=100, dt=0.001):
        """
        模拟可视化控制器
        :param simulator: IPS模拟器对象
        :param particle_system: 粒子系统对象
        :param rad: 约束半径
        :param draw_interval: 绘图间隔步数
        :param dt: 时间步长
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
        """初始化绘图元素"""
        # 位置子图
        
        self.ax1 = self.fig.add_subplot(121)
        self.scat = self.ax1.scatter([], [], s=50, edgecolors='k')
        self._draw_constraint()
        self.ax1.set_xlim(-self.rad-2, self.rad+2)
        self.ax1.set_ylim(-self.rad-2, self.rad+2)
        self.ax1.set_aspect('equal')
        self.ax1.set_title(f'Particle Positions (R={self.rad})')
        
        # 能量子图
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
        """绘制约束区域"""
        constraint = plt.Circle((0,0), self.rad, color='r', fill=False, 
                               linestyle='--', alpha=0.7)
        self.ax1.add_artist(constraint)
        
    def _init_energy_containers(self):
        """初始化能量数据容器"""
        self.energy_history = {
            'total': [],
            'kinetic': [],
            'potential': []
        }
        
    def add_callback(self, callback_func):
        """
        添加自定义回调函数
        :param callback_func: 函数格式 func(simulator, current_step)
        """
        self.callbacks.append(callback_func)
        
    def _update_energy_plot_limits(self):
        """动态更新能量图Y轴范围"""
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
        """计算各能量分量"""
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
        """更新可视化元素"""
        # 更新位置
        positions = self.p.get_positions()
        self.scat.set_offsets(np.c_[positions[0], positions[1]])
        
        # 更新能量曲线
        energies = self._calculate_energies()
        for key in self.energy_history:
            self.energy_history[key].append(energies[key])
            self.lines[key].set_data(
                range(len(self.energy_history[key])),
                self.energy_history[key]
            )
            
        # 自动调整坐标轴
        self._update_energy_plot_limits()
        self.ax2.set_xlim(0, len(self.energy_history['total']))
        
    def run_animation(self, total_steps):
        print(f"Running animation for {total_steps} steps")
        """运行动画"""
        def _frame_update(frame):
            # 执行自定义回调
            for cb in self.callbacks:
                cb(self.simulator, frame*self.draw_interval)
                
            self._update_plots(frame)
            # 推进模拟
            self.simulator.integrate_n_steps(self.dt, self.draw_interval)

            return self.scat, *self.lines.values()
        
        ani = animation.FuncAnimation(
            self.fig, _frame_update,
            frames=range(0, total_steps, self.draw_interval), init_func=self._init_plots,
        )
        print(f"Output number of frames: {total_steps // self.draw_interval}")
        # ani.save('test.gif', writer='ffmpeg', fps=10)
        json = ani.to_jshtml()
        plt.close()
        return HTML(json)
    
    def run_animation_to_gif(self, total_steps, gif_path, fps=10):
        """
        边模拟边写入 GIF 文件
        :param total_steps: 总模拟步数
        :param gif_path: 输出 GIF 文件路径
        :param fps: 帧率（每秒帧数）
        """
        print(f"Running animation for {total_steps} steps, saving to {gif_path}")
        
        # 创建临时目录存储单帧图像
        temp_dir = tempfile.mkdtemp()
        frame_files = []
        self._init_plots()
        try:
            # 逐帧生成并保存
            for frame_idx, step in enumerate(range(0, total_steps, self.draw_interval)):
                # 执行自定义回调
                for cb in self.callbacks:
                    cb(self.simulator, step)

                # 更新绘图
                self._update_plots(frame_idx)
                
                # 推进模拟
                self.simulator.integrate_n_steps(self.dt, self.draw_interval)

                # 保存当前帧为临时文件
                frame_path = os.path.join(temp_dir, f"frame_{frame_idx:04d}.png")
                self.fig.savefig(frame_path, dpi=100)
                frame_files.append(frame_path)

                # print(f"Frame {frame_idx} saved to {frame_path}")

            # 将所有帧合并为 GIF
            # self._create_gif(frame_files, gif_path, fps)

        finally:
            print(f"Output number of frames: {len(frame_files)}")
            # out put a gif
            self._create_gif(frame_files, gif_path, fps)
            # 清理临时文件
            for frame_file in frame_files:
                os.remove(frame_file)
            os.rmdir(temp_dir)

    def _create_gif(self, frame_files, gif_path, fps):
        """
        将帧图像合并为 GIF
        :param frame_files: 帧文件路径列表
        :param gif_path: 输出 GIF 文件路径
        :param fps: 帧率
        """
        images = []
        for frame_file in frame_files:
            img = Image.open(frame_file)
            images.append(img)

        # 保存为 GIF
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            optimize=True,
            duration=1000 // fps,  # 每帧持续时间（毫秒）
            loop=0  # 无限循环
        )
        print(f"GIF saved to {gif_path}")
