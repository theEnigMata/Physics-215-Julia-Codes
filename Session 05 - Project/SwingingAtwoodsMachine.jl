using Plots
const G = 9.80665   # Acceleration due to gravity g (in m/s^2)

"""
    get_dstate(s, μ)

Computes for change in initial state `s` in system given mass ratio `μ`.
"""
function get_dstate(s, μ)
    ṡ = copy(s)

    r = s[1]
    θ = s[2]
    vr = s[3]
    vθ = s[4]

    ṡ[1] = vr
    ṡ[2] = vθ
    ṡ[3] = (r*vθ^2.0 + G*(cos(θ) - μ)) / (1.00 + μ)
    ṡ[4] = - (2.00*vr*vθ + G*sin(θ)) / r
    return ṡ
end

"""
    RK4(motion, init_state, dt, μ)

Numerically integrates system for time interval `dt` given initial state `init_state` and mass ratio `μ`, stores everything in `motion`.
"""
function RK4!(motion, init_state, dt, μ)
    for i in 1:size(motion)[2]
        motion[:,i] = init_state

        ds1 = dt * get_dstate(init_state, μ)
        ds2 = dt * get_dstate(init_state + 0.5*ds1, μ)
        ds3 = dt * get_dstate(init_state + 0.5*ds2, μ)
        ds4 = dt * get_dstate(init_state + ds3, μ)
        init_state += (ds1 + 2*ds2 + 2*ds3 + ds4) / 6

        if init_state[1] <= 0.00 || init_state[1] >= 7.00
            break
        end
    end
end

"""
    plotter(motion, time, μ)

Plots the trajectory of the small mass from `motion` at all times `time` given mass ratio `μ`.
"""
function plotter(motion, time, μ)
    r = motion[1,:]
    θ = motion[2,:]
    x = r.*sin.(θ)
    y = - r.*cos.(θ)

    plt = plot(x, y
                ,label="μ = $(μ)"
                ,linecolor=:green
    )
    plot!(plt
        ,xaxis="x-position"
        ,yaxis="y-position"
        ,legend=:bottomright
        ,aspect_ratio=:equal
    )

    display(plt);
end

"""
    SwingingAtwoodsMachine(μ::Float64=1.1185, t_fin::Float64=12.500)

Given mass ratio `μ` and final time for integration `t_fin`, it integrates the system, and plots the trajectory of the resulting motion.
"""
function SwingingAtwoodsMachine(μ::Float64=1.1185, t_fin::Float64=12.500)
    dt = 10.0^(-3.0)
    time = 0.0:dt:t_fin
    sys_motion = zeros(4, length(time))

    init_state = [1.00, pi/2.00, 0.00, 0.00]
    RK4!(sys_motion, init_state, dt, μ)

    plotter(sys_motion, time, μ)
end

"""
    benchmark_func(μ::Float64=1.1185, t_fin::Float64=12.500)

Does the same thing as SwingingAtwoodsMachine() without plotting the trajectory of the resulting motion. Used mainly for benchmarking.
"""
function benchmark_func(μ::Float64=1.1185, t_fin::Float64=12.500)
    dt = 10.0^(-3.0)
    time = 0.0:dt:t_fin
    sys_motion = zeros(4, length(time))

    init_state = [1.00, pi/2.00, 0.00, 0.00]
    RK4!(sys_motion, init_state, dt, μ)
end