"""@author: Chris Sevilla"""

from numpy import arange, zeros, array, pi, sin, cos
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt

def get_dstate(s, mu):
    ds = zeros(len(s))
    
    r = s[0]
    theta = s[1]
    vr = s[2]
    vtheta = s[3]
    
    ds[0] = vr
    ds[1] = vtheta
    ds[2] = (r*vtheta**2 + G*(cos(theta) - mu)) / (1.00 + mu)
    ds[3] = - (2.00*vr*vtheta + G*sin(theta)) / r
    
    return ds


def RK4(motion, init_state, dt, mu):
    for i in range(motion.shape[0]):
        motion[i] = init_state
        
        ds1 = dt * get_dstate(init_state, mu)
        ds2 = dt * get_dstate(init_state + 0.5*ds1, mu)
        ds3 = dt * get_dstate(init_state + 0.5*ds2, mu)
        ds4 = dt * get_dstate(init_state + ds3, mu)
        init_state += (ds1 + 2*ds2 + 2*ds3 + ds4) / 6

        if init_state[0] <= 0.00 or init_state[0] >= 7.00:
            break
    
    return motion

def plotter(motion, time, mu):
    r = motion[:,0]
    theta = motion[:,1]
    
    x = r*sin(theta)
    y = -r*cos(theta)
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel(r"$x$-position")
    ax.set_ylabel(r"$y$-position")
    ax.plot(x, y)
    plt.tight_layout()
    plt.show()


def SwingingAtwoodsMachine(mu = 1.1185, t_fin = 12.50):
    dt = 10.0**(-2.0)
    time = arange(0.0, t_fin, dt)
    sys_motion = zeros((len(time), 4))
    
    init_state = array([1.00, pi/2.00, 0.00, 0.00])
    sys_motion = RK4(sys_motion, init_state, dt, mu)
    
    plotter(sys_motion, time, mu)

def benchmark_func(mu = 1.1185, t_fin = 12.50):
    dt = 10.0**(-2.0)
    time = arange(0.0, t_fin, dt)
    sys_motion = zeros((len(time), 4))
    
    init_state = array([1.00, pi/2.00, 0.00, 0.00])
    sys_motion = RK4(sys_motion, init_state, dt, mu)

G = 9.80665
def main():
    plt.close("all")

    # For running the program
    SwingingAtwoodsMachine()

    # For benchmarking the program
    import timeit
    num_of_iterations = 766
    total_time = timeit.Timer(benchmark_func).timeit(number=num_of_iterations)
    print("Average time per iteration is"+str(total_time/num_of_iterations)+" s.")
    
if __name__ == '__main__':
    main()